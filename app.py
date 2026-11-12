# app.py
import os
import traceback
import numpy as np
import pandas as pd

# 无头环境绘图（服务器/宝塔/容器下）
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.optimize import brentq
import gradio as gr


# Corrected B-atom mass balance; see README.md.
from boron_speciation_core import solve_boric_acid, molecular_concentrations, BORON_NUMBERS

def total_boron_from_x_y(x, y):
    return float(BORON_NUMBERS @ molecular_concentrations(x, y))

def solve_x_for_y(C_total, y):
    return solve_boric_acid(C_total, -np.log10(y))

def species_from_xy(x, y, C_total):
    fractions = BORON_NUMBERS * molecular_concentrations(x, y) / C_total
    k = fractions[[5, 4, 3, 2, 1]]
    return (*k, *np.cumsum(k))

def parse_concentrations(text: str):
    """支持：'0.02,0.05,0.08' 或 '0.01:0.10:0.01' 或 '0.08'。"""
    t = text.strip()
    if ":" in t:
        s, e, st = map(float, t.split(":"))
        if st <= 0:
            raise ValueError("切片步长必须 > 0")
        n = int(round((e - s) / st)) + 1
        vals = [round(s + i * st, 12) for i in range(max(n, 0))]
    elif "," in t:
        vals = [float(v) for v in t.split(",") if v.strip()]
    else:
        vals = [float(t)]
    vals = sorted(set(v for v in vals if v > 0))
    if not vals:
        raise ValueError("没有有效浓度")
    return vals


# ---------- Gradio 回调 ----------
def compute_and_export(conc_text, ph_min, ph_max, ph_points, show_plot, progress=gr.Progress(track_tqdm=True)):
    try:
        concentrations = parse_concentrations(conc_text)
        if ph_min >= ph_max:
            raise ValueError("pH 最小值必须小于最大值")
        ph_points = int(ph_points)
        if ph_points < 5:
            raise ValueError("pH 采样点数太少（建议 ≥ 50）")

        # 数值计算/积分用的网格（仅在设定区间）
        pH_grid = np.linspace(ph_min, ph_max, ph_points)
        y_grid  = 10.0 ** (-pH_grid)

        summary_rows = []
        # 写到工作目录，便于下载
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
                    "B(OH)3_fraction": xs / C,
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

        # 绘图：范围固定 pH=2~14；用虚线标注积分区间；标题仅“Cumulative Y”
        fig = None
        if show_plot:
            pH_plot_min, pH_plot_max = 2.0, 14.0
            pH_plot = np.linspace(pH_plot_min, pH_plot_max, 600)
            y_plot  = 10.0 ** (-pH_plot)

            C_last = concentrations[-1]
            Yplot = []
            for yy in y_plot:
                x = solve_x_for_y(C_last, yy)
                *_, Y1, Y2, Y3, Y4, Y5 = species_from_xy(x, yy, C_last)
                Yplot.append([Y1, Y2, Y3, Y4, Y5])
            Yplot = np.array(Yplot)

            fig = plt.figure(figsize=(7, 5))
            labels = ["Y1", "Y2", "Y3", "Y4", "Y5"]
            for i, lab in enumerate(labels):
                plt.plot(pH_plot, Yplot[:, i], label=lab, linewidth=1.2)

            plt.axvline(x=ph_min, linestyle="--")
            plt.axvline(x=ph_max, linestyle="--")
            plt.xlim(pH_plot_min, pH_plot_max)
            plt.title("Cumulative Y")
            plt.xlabel("pH")
            plt.ylabel("Y (cumulative)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()

        return tmp_xlsx, pd.DataFrame(summary_rows), fig

    except Exception as e:
        # 打印到日志并给用户提示
        print("ERROR in compute_and_export:\n", traceback.format_exc(), flush=True)
        gr.Warning(f"运行出错：{e}")
        # 返回空结果占位，避免前端组件错误
        empty_xlsx = None
        empty_table = pd.DataFrame({"message": ["error"], "detail": [str(e)]})
        empty_plot = plt.figure(figsize=(4, 2))
        plt.text(0.5, 0.5, "Error", ha="center"); plt.axis("off")
        return empty_xlsx, empty_table, empty_plot


# ---------- Gradio UI ----------
with gr.Blocks(
    title="硼的形式计算器（多浓度批量）",
    theme=gr.themes.Soft(primary_hue="blue", secondary_hue="cyan", neutral_hue="gray")
) as demo:
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
    # 本地/服务器直接运行
    demo.launch(
        server_name="0.0.0.0",
        server_port=int(os.getenv("PORT", 7860)),
        share=False,
        show_error=True
    )
