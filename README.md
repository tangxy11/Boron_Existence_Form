# 🧪 硼的存在形式计算器（Boron Existence Form Calculator）

**作者（Author）:** 汤昕煜 (Tang Xinyu)  
**单位（Affiliation）:** Qinghai University  
**License:** MIT License  

---

## 项目简介（Project Overview）

本项目是一个基于 **Python + Gradio** 构建的科研辅助工具，用于计算不同浓度条件下硼元素在水溶液体系中的存在形式（包括多种阴离子物种的比例分布），并在指定的 pH 区间内进行数值积分与可视化分析。

本工具可自动批量计算多组浓度下的平衡组成，输出结果包括：
- 各物种分数（k₁ ~ k₅）
- 累积分布函数（Y₁ ~ Y₅）
- 在指定 pH 区间的积分值（∫kᵢ dpH）
- 可视化图表（pH 2–14 范围下的累计分布曲线）

最终结果将导出为 Excel 文件，包含多个 sheet（每个浓度一个 sheet + 汇总表）。

---

## Project Description

This project provides a scientific computation tool built with **Python** and **Gradio**, designed to simulate the **speciation of boron** in aqueous systems under different concentrations and pH ranges.

It automatically performs numerical solutions and integrations for various boron-containing species and provides:
- Species fractions (k₁–k₅)
- Cumulative species ratios (Y₁–Y₅)
- Integrated values of ∫kᵢ dpH over the selected pH interval
- A full-range visualization (pH 2–14) with shaded integration boundaries

All calculated data are exported into an Excel workbook, with each concentration stored in a separate sheet, plus a summary table of integrated results.

---

## ⚙️ 核心特性（Key Features）

- ✅ 多浓度批量计算（支持 `0.01:0.10:0.01` 或 `0.02,0.05,0.08` 等格式）  
- ✅ 任意 pH 区间积分计算（例如 8.6–10.1）  
- ✅ 自动数值求解（使用 `scipy.optimize.brentq`）  
- ✅ 全范围可视化绘图（pH = 2–14，带积分区间虚线）  
- ✅ 结果导出为 Excel 文件（包含每个浓度与汇总）  
- ✅ 适配无图形界面的服务器环境（`matplotlib.use("Agg")`）

---

## 💻 使用方法（Usage）

### 1️⃣ 安装依赖
```bash
pip install -r requirements.txt

### 2️⃣ 运行程序
```bash
python app.py

### 3️⃣ 打开网页
```bash
https://huggingface.co/spaces/tangxinyu02/Boron_existence_form
即可进入交互式界面

---

📊 输出说明（Output Explanation）

输出结果包括：

✅Excel 文件：每个浓度对应一个 sheet，包含 pH、k₁–k₅、Y₁–Y₅ 等数据；另有汇总页（integrals）。

✅积分表格：∫kᵢ dpH 的结果表。

✅绘图：展示 pH = 2–14 范围内的累计分布（Cumulative Y）曲线，灰色虚线标注积分区间。

---

依赖环境（Dependencies）

| Package    | Version | Description |
| ---------- | ------- | ----------- |
| numpy      | 1.24.4  | 数值计算        |
| scipy      | 1.10.1  | 非线性方程求解     |
| pandas     | 2.0.3   | 数据处理与导出     |
| matplotlib | 3.7.2   | 绘图          |
| openpyxl   | latest  | Excel 导出支持  |
| gradio     | ≥4.0.0  | Web 界面构建    |

---

学术用途（For Research Use）

本项目旨在为科研人员提供一种快速、可视化的硼物种分布计算方法。

算法基于经典的酸碱平衡模型，适用于溶液化学、环境化学、材料吸附机理等研究。
