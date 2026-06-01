# Extract Peptide — 纳米孔阻断信号抽取流水线

基于大语言模型/视觉模型（豆包 / 火山方舟）从纳米孔单分子检测论文中自动抽取阻断信号关系的 Python 工具链。

## 概述

本工具链读取学术论文（PDF、DOCX、TXT），从中抽取以下结构化信息：

- **待测物（analyte）**：肽段、蛋白质、化合物等
- **纳米孔（nanopore）**：α-溶血素、aerolysin、MspA、ClyA、固态纳米孔等，含突变体/变体
- **阻断电流指标**：Ib/I0、I/I0、残余电流（residual current）、相对电流、阻断幅度、ΔI
- **阻断时间指标**：tD、log(tD)、停留时间（dwell time）、驻留时间（residence time）、易位时间（translocation time）
- **实验条件**：电压、缓冲液、电解质、pH、温度、加样侧（cis/trans）、待测物浓度
- **图表来源**：每个数据点对应的图号（figure）和子图（panel）

## 脚本说明

| 脚本 | 功能 |
|------|------|
| `batch_extract_nanopore_signal_relations.py` | **批量抽取主脚本。** 读取各子文件夹中的论文，对文本进行分块、评分，调用 LLM 抽取阻断信号关系，输出每篇论文的 Excel 文件及批量汇总。 |
| `filter_signal_data.py` | **后处理过滤。** 按目标浓度（5/10/15/20 μM）和电压（20 mV 的整数倍）过滤数据，清理 figure↔panel 对应关系。 |
| `read_panel_data.py` | **VLM 面板读数。** 使用视觉模型在 PDF 页面中定位子图、裁剪、评估哪个面板最适合读取 dwell_current/dwell_time，并直接从图中提取原始数值。 |

## 整体架构

```
paper1/                          ← 每篇论文一个文件夹
├── paper1.pdf                   ← 论文正文 PDF
├── paper1_si.pdf                ← 补充材料 PDF
├── paper1_nanopore_blockade_signal_relations.xlsx        ← LLM 抽取结果
├── paper1_nanopore_blockade_signal_relations_filtered.xlsx ← 过滤后结果
└── paper1_nanopore_blockade_signal_relations_filtered_vlm_evaluated.xlsx ← VLM 评估后结果
```

### 数据流

```
PDF / DOCX / TXT 论文
       │
       ▼
batch_extract_nanopore_signal_relations.py   ← LLM 文本抽取 + 可选 VLM 页面识别
       │
       ▼
*_nanopore_blockade_signal_relations.xlsx    ← 原始抽取记录
       │
       ▼
filter_signal_data.py                        ← 按浓度和电压过滤
       │
       ▼
*_filtered.xlsx                              ← 过滤后记录
       │
       ▼
read_panel_data.py                           ← VLM 子图定位与数值读取
       │
       ▼
*_vlm_evaluated.xlsx                         ← 最终输出（含 VLM 读取的数值）
```

## 环境依赖

```bash
pip install pymupdf pandas openpyxl requests
# 可选：读取 DOCX 文件需要
pip install python-docx
```

## 环境变量

| 变量名 | 说明 |
|--------|------|
| `ARK_API_KEY` / `DOUBAO_API_KEY` / `VOLCENGINE_API_KEY` | 豆包/火山方舟 API Key，三者任选其一 |
| `DOUBAO_MODEL` / `ARK_MODEL` / `VOLCENGINE_MODEL` | 文本模型名或推理接入点 ID |
| `DOUBAO_VL_MODEL` / `ARK_VL_MODEL` | 视觉模型名或视觉推理接入点 ID（仅 `--use-vlm-pages` 时使用） |
| `ARK_BASE_URL` | API 接口地址（默认：`https://ark.cn-beijing.volces.com/api/v3/chat/completions`） |

## 使用方法

### 1. 批量抽取阻断信号关系

```bash
# 基础抽取（仅文本）
python batch_extract_nanopore_signal_relations.py "D:/Code/Python/extract_peptide"

# 启用 VLM 读取图表页面
python batch_extract_nanopore_signal_relations.py "D:/Code/Python/extract_peptide" --use-vlm-pages

# 限制每个文件夹 VLM 读取的页数
python batch_extract_nanopore_signal_relations.py "D:/Code/Python/extract_peptide" --use-vlm-pages --vlm-max-pages-per-folder 20
```

### 2. 过滤抽取数据

```bash
# 过滤所有论文文件夹
python filter_signal_data.py "D:/Code/Python/extract_peptide"

# 过滤单个 Excel 文件
python filter_signal_data.py --excel "D:/path/to/file.xlsx"
```

### 3. VLM 读取面板数据

```bash
# 处理所有论文文件夹
python read_panel_data.py "D:/Code/Python/extract_peptide"

# 处理指定文件夹并使用自定义 VLM 模型
python read_panel_data.py "D:/Code/Python/extract_peptide" --folder paper1 --vlm-model ep-xxxxxxxx
```

## 输出字段说明

Excel 输出的 `signal_relation_records` 工作表包含以下关键字段：

| 字段 | 说明 |
|------|------|
| `analyte` | 待测物名称（肽段/蛋白质/化合物） |
| `nanopore` | 纳米孔名称（aerolysin、α-hemolysin 等） |
| `nanopore_variant` | 纳米孔突变体/变体（如 T232K/K238Q） |
| `dwell_current_metric` | 电流指标类型（Ib/I0、residual current 等） |
| `dwell_time_metric` | 时间指标类型（tD、dwell time 等） |
| `dwell_current` | 提取的电流数值 |
| `dwell_time` | 提取的时间数值 |
| `voltage_mV` | 电压（mV） |
| `analyte_concentration_value` | 待测物浓度数值 |
| `electrolyte` | 电解质（如 1 M KCl） |
| `buffer` | 缓冲液（如 10 mM Tris） |
| `pH` | pH 值 |
| `figure_id_dwell_current` | 电流数据来源图号 |
| `figure_id_dwell_time` | 时间数据来源图号 |
| `panel_dwell_current` | 电流数据对应子图 |
| `panel_dwell_time` | 时间数据对应子图 |

## 关键设计决策

- **σb / current fluctuation / 噪声不作为目标信号输出** —— 这些仅在 paired_panels、panel_relation 或 evidence 中作为辅助上下文保留，不会单独生成目标记录。
- **电压是关系键，不同电压不合并** —— 每个电压值单独占一行，电压序列如 80/100/120/140/160 mV 会被自动展开。
- **浓度按待测物区分** —— 不同待测物有不同浓度时分别记录，不会合并为笼统值。
- **Figure ID 拆分为电流和时间两列** —— `figure_id_dwell_current` 和 `figure_id_dwell_time` 独立存放，因为电流和时间数据可能来自不同的图。
- **肽段编号/突变编号不作为电压** —— 解析器主动过滤名称中的数字（如 "Aβ18-26" 中的 18 和 26、"A21G" 中的 21、"T232K/K238Q" 中的 232 和 238），避免误识别为电压值。
- **面板维度优先级** —— 当同一张图中存在多种维度的面板时，优先使用一维图（直方图、柱状图、电压依赖折线图）；二维散点图、等高线图、热力图不作为目标数据来源。

## 项目结构

```
extract_peptide/
├── README.md                                     ← 本文件
├── batch_extract_nanopore_signal_relations.py    ← 批量抽取主脚本
├── filter_signal_data.py                         ← 数据过滤脚本
├── read_panel_data.py                            ← VLM 面板读取脚本
├── batch_nanopore_blockade_signal_summary.xlsx   ← 批量汇总输出
└── paper1/                                       ← 示例论文文件夹
    ├── paper1.pdf
    ├── paper1_si.pdf
    ├── paper1_nanopore_blockade_signal_relations.xlsx
    ├── paper1_nanopore_blockade_signal_relations_filtered.xlsx
    └── _panel_read_debug/                        ← 面板读取调试输出
```
