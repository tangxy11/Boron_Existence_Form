# 🧪 硼的存在形式计算器（Boron Existence Form Calculator）

**作者（Author）:** 汤昕煜 (Tang Xinyu)  
**单位（Affiliation）:** Qinghai University  
**License:** MIT License  

---

## 🇨🇳 项目简介（Project Overview）

本项目基于 **Python + Gradio** 构建，旨在计算不同浓度条件下硼元素在水溶液体系中的存在形式，包括多种阴离子物种（如 B(OH)₃、B(OH)₄⁻ 等）的分布比例。  
同时可对选定 pH 区间（如 8.6–10.1）进行数值积分，并提供可视化分析。

### ✨ 功能摘要

- 🔹 计算各物种分数（k₁ ~ k₅）  
- 🔹 计算累积分布函数（Y₁ ~ Y₅）  
- 🔹 在自定义 pH 区间计算积分值（∫kᵢ dpH）  
- 🔹 绘制 pH 2–14 的全范围可视化曲线（带积分虚线）  
- 🔹 导出 Excel 文件（每个浓度独立 sheet + 汇总表）  

---

## 💻 使用方法（Usage）

### 1️⃣ 安装依赖
```bash
pip install -r requirements.txt
```

---

### 2️⃣ 运行程序
```bash
python app.py
```

---

### 3️⃣ 打开网页
```bash
https://huggingface.co/spaces/tangxinyu02/Boron_existence_form
即可进入交互式界面
```

---

## 📊 输出说明（Output Explanation）

- 输出结果包括：

- ✅ Excel 文件：每个浓度对应一个 sheet，包含 pH、k₁–k₅、Y₁–Y₅ 等数据；另有汇总页（integrals）。

- ✅ 积分表格：∫kᵢ dpH 的结果表。

- ✅ 绘图：展示 pH = 2–14 范围内的累计分布（Cumulative Y）曲线，灰色虚线标注积分区间。

---

## 依赖环境（Dependencies）

| Package    | Version |
| ---------- | ------- |
| numpy      | 1.24.4  |
| scipy      | 1.10.1  |
| pandas     | 2.0.3   |
| matplotlib | 3.7.2   |
| openpyxl   | latest  |
| gradio     | ≥4.0.0  |

---

## 学术用途（For Research Use）

本项目旨在为科研人员提供一种快速、可视化的硼物种分布计算方法。
算法基于经典的酸碱平衡模型，适用于溶液化学、环境化学、材料吸附机理等研究。