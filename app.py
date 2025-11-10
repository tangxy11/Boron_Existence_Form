# app.py
import os, io, traceback
import numpy as np
import pandas as pd

# 无头环境绘图
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.optimize import brentq
import gradio as gr

# ---------- 核心数值模型 ----------
def total_boron_from_x_y(x, y):
    if y <= 0:
        return np.inf
    term1 = x
    term2 = (10**(-9.2))  * x / y
    term3 = (10**(-7.29)) * (x**3) / y
    term4 = (10**(-6.77)) * (x**5) / y
    term5 = (10**(-14.5)) * (x**4) / (y**2)
    term6 = (10**(-16.3)) * (x**3) / (y**2)
    return term1 + term2 + term3 + term4 + term5 + term6

def solve_x_for_y(C_total, y):
    a = 1e-16
    b = max(C_total, 1e-8)
    fa = total_boron_from_x_y(a, y) - C_total
    fb = total_boron_from_x_y(b, y) - C_total
    expand = 0
    while fa * fb > 0 and expand < 40:
        b *= 2.0
        fb = total_boron_from_x_y(b, y) - C_total
        expand += 1
    if fa * fb > 0:
        return min(b, C_total)
    root = brentq(lambda xx: total_boron_from_x_y(xx, y) - C_total, a, b, maxiter=200, xtol=1e-12)
    return root

def species_from_xy(x, y, C_total):
    a = (10**(-9.2))  * (x / y)
    b = (10**(-7.29)) * (x**3 / y)
    c = (10**(-6.77)) * (x**5 / y)
    d = (10**(-14.5)) * (x**4 / (y**2))
    e = (10**(-16.3)) * (x**3 / (y**2))

    k1 = e / C_total
    k2 = d / C_total
    k3 = c / C_total
    k4 = b / C_total
    k5 = a / C_total

    Y1 = k1
    Y2 = k1 + k2
    Y3 = k1 + k2 + k3
    Y4 = k1 + k2 + k3 + k4
    Y5 = k1 + k2 + k3 + k4 + k5
    return k1, k2, k3, k4, k5, Y1, Y2, Y3, Y4, Y5

def parse_concentrations(text):
    text = text.strip()
    if ":" in text:
        s, e, st = map(float, text.split(":"))
        if st <= 0:
            raise ValueError("切片步长必须 > 0")
        n = int(round((e - s) / st)) + 1
        vals = [round(s + i * st, 12) for i in range(max(n, 0))]
    elif "," in text:
        vals = [float(v) for v in text.split(",") if v.strip()]
    else:
        vals = [float(text)]
    vals = sorted(set(v for v in vals if v > 0))
    if not vals:
        raise ValueError("没有有效浓度")
    return vals

def compute_and_export(conc_text, ph_min, ph_max, ph_points, show_plot, progress=gr.Progress(track_tqdm=True)):
    try:
        concentrations = parse_concentrations(conc_text)
        if ph_min >= ph_max:
            raise ValueError("pH 最小值必须小于最大值")
        ph_points = int(ph_points)
        if ph_points < 5:
            raise ValueError("pH 采样点数太少（建议 ≥ 50）")

        pH_grid = np.linspace(ph_min, ph_max, ph_points)
        y_grid  = 10.0 ** (-pH_grid)

        summary_rows = []
        # 写到临时文件（某些 gradio 版本不接受 bytes tuple）
        tmp_xlsx = os.path.join(os.getcwd(), "boron_forms.xlsx")

        with pd.ExcelWriter(tmp_xlsx, engine="openpyxl") as writer:
            for i, C in enumerate(concentrations, 1):
                progress((i-1)/len(concentrations), desc=f"正在计算浓度 {C:g}")
                xs, ks_arr, Ys_arr = [], [], []
                for y in y_grid:
                    x = solve_x_for_y(C, y)
                    xs.append(x)
                    ks = species_from_xy(x, y, C)
                    ks_arr.append(ks[:5])
                    Ys_arr.append(ks[5:])

                xs = np.array(xs)
                ks_arr = np.array(ks_arr)   # (N,5)
                Ys_arr = np.array(Ys_arr)   # (N,5)

                integrals = np.trapz(ks_arr, pH_grid, axis=0)
                summary_rows.append({
                    "Concentration": C,
                    "∫k1 dpH": float(integrals[0]),
                    "∫k2 dpH": float(integrals[1]),
                    "∫k3 dpH": float(integrals[2]),
                    "∫k4 dpH": float(integrals[3]),
                    "∫k5 dpH": float(integrals[4]),
                    "pH_min": ph_min, "pH_max": ph_max, "points": ph_points
                })

                df = pd.DataFrame({
                    "pH": pH_grid,
                    "x_solved": xs,
                    "k1": ks_arr[:,0],
                    "k2": ks_arr[:,1],
                    "k3": ks_arr[:,2],
                    "k4": ks_arr[:,3],
                    "k5": ks_arr[:,4],
                    "Y1": Ys_arr[:,0],
                    "Y2": Ys_arr[:,1],
                    "Y3": Ys_arr[:,2],
                    "Y4": Ys_arr[:,3],
                    "Y5": Ys_arr[:,4],
                })
                sheet_name = f"C={C:g}"[:31]
                df.to_excel(writer, index=False, sheet_name=sheet_name)

            pd.DataFrame(summary_rows).to_excel(writer, index=False, sheet_name="integrals")

        # 绘图（可选）
                # 绘图（可选）：曲线覆盖 pH = 2 ~ 14，积分区间用虚线标示
        fig = None
        if show_plot:
            # 1) 用宽范围的 pH 网格只用于绘图
            pH_plot_min, pH_plot_max = 2.0, 14.0
            pH_plot = np.linspace(pH_plot_min, pH_plot_max, 600)   # 可调，比如 400/800
            y_plot  = 10.0 ** (-pH_plot)

            # 2) 画最后一个浓度的累计曲线
            C = concentrations[-1]
            Yplot = []
            for yy in y_plot:
                x = solve_x_for_y(C, yy)
                *_, Y1, Y2, Y3, Y4, Y5 = species_from_xy(x, yy, C)
                Yplot.append([Y1, Y2, Y3, Y4, Y5])
            Yplot = np.array(Yplot)  # shape: (N,5)

            # 3) 画图 & 标出积分区间虚线
            fig = plt.figure(figsize=(7,5))
            labels = ["Y1","Y2","Y3","Y4","Y5"]
            for i, lab in enumerate(labels):
                plt.plot(pH_plot, Yplot[:, i], label=lab, linewidth=1.2)

            # 标示积分区间
            plt.axvline(x=ph_min, linestyle="--")
            plt.axvline(x=ph_max, linestyle="--")

            # 边界与样式
            plt.xlim(pH_plot_min, pH_plot_max)
            plt.title(f"Cumulative Y (Integration interval {ph_min}–{ph_max}, C={C:g})")
            plt.xlabel("pH"); plt.ylabel("Y (cumulative)")
            plt.grid(True); plt.legend(); plt.tight_layout()


        return tmp_xlsx, pd.DataFrame(summary_rows), fig

    except Exception as e:
        # 打印到 Runtime logs，界面弹警告
        print("ERROR in compute_and_export:\n", traceback.format_exc(), flush=True)
        gr.Warning(f"运行出错：{e}")
        # 返回三个“空”输出，避免 Gradio 输出组件报红
        empty_xlsx = None
        empty_table = pd.DataFrame({"message": ["error"], "detail": [str(e)]})
        empty_plot = plt.figure(figsize=(4,2)); plt.text(0.5,0.5,"Error", ha="center"); plt.axis("off")
        return empty_xlsx, empty_table, empty_plot

# ---------- Gradio UI ----------
with gr.Blocks(title="硼的形式计算器（多浓度批量）") as demo:
    gr.Markdown("## 硼的存在形式计算（多浓度批量）")
    with gr.Row():
        conc = gr.Textbox(label="浓度（示例：0.02,0.05,0.08 或 0.01:0.10:0.01）", value="0.08")
    with gr.Row():
        ph_min = gr.Number(label="pH 最小值", value=8.6)
        ph_max = gr.Number(label="pH 最大值", value=10.1)
        ph_pts = gr.Number(label="pH 采样点数（建议≥300）", value=300, precision=0)
        show_plot = gr.Checkbox(label="显示示意图（最后一个浓度）", value=True)

    run_btn = gr.Button("运行计算")
    excel_out = gr.File(label="下载 Excel（多 sheet + 汇总）", file_count="single")
    table_out = gr.Dataframe(label="积分汇总（∫k_i dpH）", wrap=True)
    plot_out = gr.Plot(label="累计曲线示意（Y1..Y5）")

    run_btn.click(
        fn=compute_and_export,
        inputs=[conc, ph_min, ph_max, ph_pts, show_plot],
        outputs=[excel_out, table_out, plot_out]
    )

if __name__ == "__main__":
    demo.launch()
