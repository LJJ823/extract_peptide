#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
read_panel_data.py

读取过滤后的 Excel，按 (analyte, nanopore_variant, electrolyte, pH,
analyte_concentration_value, voltage_raw) 分组后，对每组内的所有候选
(figure_id_dwell_current, panel_dwell_current) 面板逐一评估，选最优面板提取
dwell_current / dwell_time。

用法：
    python read_panel_data.py "D:/Code/Python/extract_peptide/paper1"
    python read_panel_data.py "D:/Code/Python/extract_peptide" --folder paper1
    python read_panel_data.py "D:/Code/Python/extract_peptide" --folder paper1 --vlm-model ep-xxx
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests

try:
    import fitz  # PyMuPDF
except Exception:
    raise RuntimeError("请先安装 PyMuPDF：pip install pymupdf")

# =============================================================================
# 1. 配置
# =============================================================================

DEFAULT_BASE_URL = "https://ark.cn-beijing.volces.com/api/v3/chat/completions"
SKIP_FOLDERS = {"__pycache__", ".git", ".idea", ".vscode", ".ipynb_checkpoints"}
FIGURE_ID_RE = re.compile(r"(?:Fig\.?|Figure)\s*(S?\d+[A-Za-z]?)", re.I)

GROUP_KEY_COLS = [
    "analyte", "nanopore_variant", "electrolyte", "pH",
    "analyte_concentration_value", "voltage_raw",
]


# =============================================================================
# 2. VLM 客户端
# =============================================================================

class DoubaoArkVLM:
    def __init__(
        self,
        api_key: Optional[str] = None,
        model: Optional[str] = None,
        base_url: Optional[str] = None,
        timeout: int = 180,
        max_retries: int = 4,
    ) -> None:
        self.api_key = (
            api_key
            or os.getenv("ARK_API_KEY", "").strip()
            or os.getenv("DOUBAO_API_KEY", "").strip()
            or os.getenv("VOLCENGINE_API_KEY", "").strip()
        )
        self.model = (
            model
            or os.getenv("DOUBAO_VL_MODEL", "").strip()
            or os.getenv("ARK_VL_MODEL", "").strip()
            or os.getenv("DOUBAO_MODEL", "").strip()
            or os.getenv("ARK_MODEL", "").strip()
            or "doubao-seed-2-0-pro-260215"
        )
        self.base_url = (
            base_url
            or os.getenv("ARK_BASE_URL", "").strip()
            or os.getenv("DOUBAO_BASE_URL", "").strip()
            or DEFAULT_BASE_URL
        )
        self.timeout = timeout
        self.max_retries = max_retries
        if not self.api_key:
            raise ValueError(
                "未检测到 API Key。请设置 ARK_API_KEY / DOUBAO_API_KEY"
            )

    def chat(
        self, messages: List[Dict[str, Any]],
        temperature: float = 0.05, max_tokens: int = 4000,
    ) -> str:
        payload = {
            "model": self.model, "messages": messages,
            "temperature": temperature, "max_tokens": max_tokens,
        }
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }
        last_err: Optional[Exception] = None
        for attempt in range(1, self.max_retries + 1):
            try:
                resp = requests.post(
                    self.base_url, headers=headers, json=payload,
                    timeout=self.timeout,
                )
                if resp.status_code >= 400:
                    raise RuntimeError(f"HTTP {resp.status_code}: {resp.text[:1200]}")
                return resp.json()["choices"][0]["message"]["content"]
            except Exception as exc:
                last_err = exc
                if attempt == self.max_retries:
                    break
                print(f"  [重试 {attempt}/{self.max_retries}] {exc}")
                time.sleep(min(20, 2 ** attempt))
        raise RuntimeError(f"VLM 调用失败：{last_err}")


# =============================================================================
# 3. PDF 操作
# =============================================================================

def render_page_to_base64(pdf_path: Path, page_no: int, zoom: float = 2.0) -> str:
    doc = fitz.open(pdf_path)
    try:
        page = doc.load_page(page_no - 1)
        pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        return base64.b64encode(pix.tobytes("png")).decode("ascii")
    finally:
        doc.close()


def render_region_to_base64(
    pdf_path: Path, page_no: int,
    rect: Tuple[float, float, float, float], zoom: float = 3.0,
) -> str:
    doc = fitz.open(pdf_path)
    try:
        page = doc.load_page(page_no - 1)

        # 确保 rect 有效：x0<x1, y0<y1，最小尺寸 20pt
        x0, y0, x1, y1 = float(rect[0]), float(rect[1]), float(rect[2]), float(rect[3])
        x0, x1 = min(x0, x1), max(x0, x1)
        y0, y1 = min(y0, y1), max(y0, y1)

        page_w = page.rect.width
        page_h = page.rect.height
        x0 = max(0, x0)
        y0 = max(0, y0)
        x1 = min(page_w, x1)
        y1 = min(page_h, y1)

        MIN_SIZE = 20
        if x1 - x0 < MIN_SIZE:
            cx = (x0 + x1) / 2
            x0 = max(0, cx - MIN_SIZE / 2)
            x1 = min(page_w, cx + MIN_SIZE / 2)
        if y1 - y0 < MIN_SIZE:
            cy = (y0 + y1) / 2
            y0 = max(0, cy - MIN_SIZE / 2)
            y1 = min(page_h, cy + MIN_SIZE / 2)

        pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), clip=fitz.Rect(x0, y0, x1, y1), alpha=False)
        return base64.b64encode(pix.tobytes("png")).decode("ascii")
    finally:
        doc.close()


def get_page_size(pdf_path: Path, page_no: int) -> Tuple[float, float]:
    doc = fitz.open(pdf_path)
    try:
        rect = doc.load_page(page_no - 1).rect
        return rect.width, rect.height
    finally:
        doc.close()


# =============================================================================
# 4. Figure / Panel 定位
# =============================================================================

def is_supplementary_figure(figure_id: str) -> bool:
    return bool(re.search(r"Fig(?:ure)?\.?\s*S\d+", str(figure_id), re.I))


def find_figure_page(pdf_path: Path, figure_id: str) -> Optional[int]:
    """在 PDF 中定位 figure_id 所在页码（1-based），多模式评分选最优。

    不再取第一个出现位置（那可能是目录/交叉引用），而是对每页评分，
    返回分最高的页面。
    """
    m = FIGURE_ID_RE.search(str(figure_id))
    if not m:
        return None
    fig_num = m.group(1)
    escaped = re.escape(fig_num)

    # 按匹配质量递减排列：(pattern, score)
    patterns: list[tuple[re.Pattern, int]] = [
        # 0: 标准 caption — "Figure 5." / "Fig. S6:" / "Figure 5 (a)"
        (re.compile(
            rf"(?:^|(?<=[\n.]))\s*(?:Fig\.?|Figure)\s*{escaped}\b"
            rf"(?:\s*[\.:\–\–\-]|\s+[A-Z]|\s*\()",
            re.I,
        ), 100),
        # 1: "Supplementary Figure S6" / "Supporting Figure S6"
        (re.compile(
            rf"(?:Supplementary|Supporting)\s+(?:Fig\.?|Figure)\s*{escaped}\b",
            re.I,
        ), 95),
        # 2: 行首 "Figure S6"
        (re.compile(
            rf"^\s*(?:Supplementary\s+)?(?:Fig\.?|Figure)\s*{escaped}\b",
            re.I | re.MULTILINE,
        ), 80),
    ]
    # 宽松兜底
    loose_pat = re.compile(
        rf"(?:Fig\.?|Figure)\s*{escaped}\b", re.I,
    )

    # panel 描述 → 很可能是真正的 figure 页面
    panel_hint = re.compile(r"\([a-h]\)|panel\s+[a-h]", re.I)
    # 目录/参考文献/缩略图页，降权
    toc_ref_hint = re.compile(
        r"(?:table\s+of\s+contents?|references?\s*$|bibliography|"
        r"supplementary\s+figures?\s*$|list\s+of\s+figures?)",
        re.I,
    )

    doc = fitz.open(pdf_path)
    try:
        scored: list[tuple[int, int]] = []
        for i in range(doc.page_count):
            text = doc.load_page(i).get_text("text") or ""
            if not text.strip():
                continue

            score = 0
            for pat, pts in patterns:
                if pat.search(text):
                    score = pts
                    break
            if score == 0:
                if loose_pat.search(text):
                    score = 30
            if score == 0:
                continue

            if panel_hint.search(text):
                score += 20
            if toc_ref_hint.search(text):
                score -= 50

            if score > 0:
                scored.append((i + 1, score))

        if not scored:
            return None

        scored.sort(key=lambda x: x[1], reverse=True)
        best_page, best_score = scored[0]
        if len(scored) > 1:
            print(f"      页面评分: {scored[:4]}")
        return best_page
    finally:
        doc.close()


def find_figure_in_pdfs(
    pdf_files: List[Path], figure_id: str,
) -> Tuple[Optional[Path], Optional[int]]:
    """在多个 PDF 中定位 figure，SI 图优先搜索 SI 文件。"""
    is_si = is_supplementary_figure(str(figure_id))
    si_kw = {"si", "supp", "supplement", "supplementary", "supporting", "esm", "appendix"}

    def _sort_key(p: Path) -> int:
        lower = p.name.lower()
        is_si_file = any(kw in lower for kw in si_kw)
        return 0 if (is_si == is_si_file) else 1

    for pdf_path in sorted(pdf_files, key=_sort_key):
        page_no = find_figure_page(pdf_path, str(figure_id))
        if page_no:
            return pdf_path, page_no
    return None, None


def vlm_identify_panels(
    vlm: DoubaoArkVLM, pdf_path: Path, page_no: int,
    target_panels: List[str],
) -> Dict[str, Tuple[float, float, float, float]]:
    """VLM 识别页面中各 panel 的位置，返回 {label: (x0,y0,x1,y1)}。"""
    b64 = render_page_to_base64(pdf_path, page_no, zoom=2.0)
    panels_str = ", ".join(target_panels)

    prompt = f"""识别这张论文页面截图中所有子图（panel）的位置。
需要找的 panel：{panels_str}

返回 JSON（百分比坐标）：
{{"panels": [
  {{"label": "a", "x0_pct": 0.05, "y0_pct": 0.1, "x1_pct": 0.45, "y1_pct": 0.5}}
]}}
只返回 JSON，不要其他文字。"""

    messages = [
        {"role": "user", "content": [
            {"type": "text", "text": prompt},
            {"type": "image_url", "image_url": {"url": f"data:image/png;base64,{b64}"}},
        ]}
    ]

    try:
        raw = vlm.chat(messages, temperature=0.05, max_tokens=2000)
    except Exception as exc:
        print(f"    VLM 识别 panel 位置失败: {exc}")
        return {}

    try:
        m = re.search(r'\{.*\}', raw, re.DOTALL)
        if not m:
            return {}
        data = json.loads(m.group())
    except Exception:
        return {}

    page_w, page_h = get_page_size(pdf_path, page_no)
    result = {}
    for info in data.get("panels", []):
        label = str(info.get("label", "")).strip().lower()
        if label not in [p.lower() for p in target_panels]:
            continue
        margin = 10
        x0 = max(0, info.get("x0_pct", 0) * page_w - margin)
        y0 = max(0, info.get("y0_pct", 0) * page_h - margin)
        x1 = min(page_w, info.get("x1_pct", 1) * page_w + margin)
        y1 = min(page_h, info.get("y1_pct", 1) * page_h + margin)
        result[label] = (x0, y0, x1, y1)
    return result


# =============================================================================
# 5. 面板评估：判断是否适合读取 dwell_current
# =============================================================================

def _is_blank_value(value: Any) -> bool:
    """判断单元格值是否为空。"""
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except Exception:
        pass
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "null", "na", "n/a"}


def _shorten_text(text: Any, max_len: int = 300) -> str:
    """压缩空白并截断，避免 prompt 过长。"""
    text = re.sub(r"\s+", " ", str(text)).strip()
    if len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


def build_group_context(
    group_df: pd.DataFrame,
    group_info: Dict[str, Any],
    max_chars: int = 2200,
) -> str:
    """把当前分组的关键信息压缩成 VLM 可用的上下文。

    只用于提示 VLM 该读哪一个 analyte / nanopore / voltage / concentration，
    不做单位标准化。
    """
    lines: List[str] = []

    main_parts = []
    for key in GROUP_KEY_COLS:
        value = group_info.get(key)
        if not _is_blank_value(value):
            main_parts.append(f"{key}={_shorten_text(value, 120)}")
    if main_parts:
        lines.append("分组条件：" + "; ".join(main_parts))

    # 优先给 VLM 看这些类型的列：条件、图号、panel、原始证据/图注等。
    include_keywords = (
        "analyte", "nanopore", "variant", "pore", "electrolyte", "ph",
        "concentration", "voltage", "figure", "fig", "panel",
        "condition", "buffer", "solution", "caption", "evidence",
        "sentence", "context", "note", "raw", "peptide", "protein",
    )
    exclude_keywords = (
        "candidate_evaluation_summary", "vlm_", "best_", "score",
    )

    for col in group_df.columns:
        col_lower = str(col).lower()
        if any(k in col_lower for k in exclude_keywords):
            continue
        if not any(k in col_lower for k in include_keywords):
            continue

        values: List[str] = []
        seen = set()
        for value in group_df[col].tolist():
            if _is_blank_value(value):
                continue
            text = _shorten_text(value, 260)
            if text in seen:
                continue
            seen.add(text)
            values.append(text)
            if len(values) >= 3:
                break

        if values:
            lines.append(f"{col}: " + " | ".join(values))

    context = "\n".join(lines)
    if len(context) > max_chars:
        context = context[: max_chars - 3] + "..."
    return context


def extract_figure_text_context(
    pdf_path: Path,
    page_no: int,
    figure_id: str,
    max_chars: int = 1600,
) -> str:
    """从 figure 所在页抽取附近文字/图注，作为读图辅助上下文。"""
    try:
        doc = fitz.open(pdf_path)
        try:
            text = doc.load_page(page_no - 1).get_text("text") or ""
        finally:
            doc.close()
    except Exception:
        return ""

    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return ""

    fig_num = _extract_figure_number(str(figure_id))
    if fig_num:
        pat = re.compile(rf"(?:Fig\.?|Figure)\s*{re.escape(fig_num)}\b", re.I)
        m = pat.search(text)
        if m:
            start = max(0, m.start() - 250)
            end = min(len(text), m.start() + max_chars)
            return _shorten_text(text[start:end], max_chars)

    return _shorten_text(text, max_chars)


def normalize_metric_type(metric_type: Any, x_label: Any = "", y_label: Any = "") -> str:
    """把 VLM 的自由输出归一到 dwell_current / dwell_time / noise / other。

    只归一化类型，不改变原始读数和单位。
    """
    text = " ".join([str(metric_type or ""), str(x_label or ""), str(y_label or "")]).lower()
    text = text.replace("−", "-").replace("–", "-").replace("—", "-")
    text = re.sub(r"\s+", " ", text)

    noise_keys = ["noise", "sigma", "σ", "rms", "fluctuation", "std", "standard deviation"]
    time_keys = [
        "dwell time", "dwell-time", "residence time", "duration",
        "translocation time", "event time", "log(td", "log t",
        "log10", "ln(td", "t_d", "td", "τ", "tau",
    ]
    current_keys = [
        "dwell_current", "current", "residual", "blockade", "blockage",
        "i/i0", "i / i0", "ib/i0", "ib / i0", "i_b/i_0",
        "i0", "i_b", "delta i", "Δi", "di/i", "conductance",
    ]

    # current fluctuation / sigma 这类优先归为 noise，避免因为含 current 被误判。
    if any(k in text for k in noise_keys):
        if not any(k in text for k in ["blockade", "blockage", "residual", "i/i0", "ib/i0", "i / i0"]):
            return "noise"
    if any(k in text for k in time_keys):
        return "dwell_time"
    if any(k in text for k in current_keys):
        return "dwell_current"
    if "noise" in text or "other" in text:
        return "noise" if "noise" in text else "other"
    return "other"


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        if _is_blank_value(value):
            return default
        return int(float(value))
    except Exception:
        return default


def evaluate_panel(
    vlm: DoubaoArkVLM,
    pdf_path: Path, page_no: int,
    panel_rect: Tuple[float, float, float, float],
    panel_label: str, figure_id: str,
    group_info: Dict[str, Any],
    target_metric: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """评估一个面板是否适合读取 dwell_current / dwell_time。

    返回原始读数和原始单位；这里不做单位统一。
    """
    b64 = render_region_to_base64(pdf_path, page_no, panel_rect, zoom=3.0)

    analyte = group_info.get("analyte", "unknown")
    voltage_raw = group_info.get("voltage_raw", "unknown")
    group_context = str(group_info.get("_context_summary", "")).strip()
    figure_text_context = extract_figure_text_context(pdf_path, page_no, figure_id)
    expected = target_metric or "dwell_current/dwell_time"

    prompt = f"""分析这张纳米孔实验图表面板，并读取最基础的原始数值。

当前任务：正在为 {expected} 选择候选 panel。
基础上下文：待测物={analyte}, Figure={figure_id}, Panel={panel_label}, 目标电压={voltage_raw}

当前分组的完整上下文：
{group_context if group_context else "无额外上下文"}

Figure 所在页文字/图注片段：
{figure_text_context if figure_text_context else "无可用 PDF 文本"}

请按以下规则判断：
1. metric_type 必须只能填下面四类之一：
   - "dwell_current": 包括 Ib/I0、I/I0、Ires/I0、residual current、remaining current、blockade current、blockage、ΔI/I0、current level 等所有电流/阻断电流/残余电流指标。
   - "dwell_time": 包括 tD、t_d、τ、dwell time、residence time、duration、translocation time、log(tD)、log10(tD) 等所有时间指标。
   - "noise": 包括 σb、sigma、RMS noise、current fluctuation 等噪声指标。
   - "other": 非上述目标指标。

2. 如果同一 panel 有多条曲线、多个柱子或多个点：
   - 优先选择与当前分组上下文最匹配的 analyte / nanopore_variant / electrolyte / pH / concentration / voltage。
   - 如果是电压依赖图，读取目标电压={voltage_raw} 对应的原始 y 值。
   - 如果是直方图/拟合曲线，读取峰位置或拟合峰值对应的原始 x 值。
   - 如果无法确定是哪一条曲线/柱子，在 value_note 中说明，不要强行编造。

3. 只读取图中的原始值和原始单位，不要换算单位，不要把百分数转小数。

4. 可读性评分 score：
   5=极清晰可较精确读取；4=较清晰可读大致值；3=能辨识但不确定；2=模糊；1=无法读取。

返回 JSON（只返回 JSON，不要解释文字）：
{{"metric_type":"dwell_current","score":4,"chart_type":"...",
 "x_axis_label":"...","y_axis_label":"...",
 "selected_series_or_bar":"...",
 "dwell_current_value":0.45,"dwell_current_unit":"dimensionless",
 "dwell_time_value":2.1,"dwell_time_unit":"log10(ms)",
 "value_note":null,"reasoning":"..."}}"""

    messages = [
        {"role": "user", "content": [
            {"type": "text", "text": prompt},
            {"type": "image_url", "image_url": {"url": f"data:image/png;base64,{b64}"}},
        ]}
    ]

    try:
        raw = vlm.chat(messages, temperature=0.05, max_tokens=2500)
    except Exception as exc:
        print(f"      VLM 评估失败: {exc}")
        return None

    try:
        m = re.search(r'\{.*\}', raw, re.DOTALL)
        if not m:
            print(f"      VLM 输出无 JSON")
            return None
        json_str = m.group()
        data = json.loads(json_str)
    except json.JSONDecodeError:
        # LLM 可能输出包含无效转义（如 \s）的 JSON，预处理修复
        try:
            fixed = re.sub(r'\\(?!["\\/bfnrtu])', r'\\\\', json_str)
            data = json.loads(fixed)
        except json.JSONDecodeError as exc2:
            print(f"      JSON 解析失败: {exc2}")
            return None
    except Exception as exc:
        print(f"      JSON 解析失败: {exc}")
        return None

    return data


# =============================================================================
# 6. Panel 名称解析
# =============================================================================

_FIG_NUM_RE = re.compile(r"(?:Figure|Fig\.?)\s*(S?\d+)", re.I)
_PANEL_ENTRY_RE = re.compile(r"^(s?\d+)([a-z])$", re.I)


def _extract_figure_number(figure_id: str) -> Optional[str]:
    """'Figure 5' -> '5', 'Figure S10' -> 'S10'"""
    m = _FIG_NUM_RE.search(str(figure_id))
    return m.group(1) if m else None


def _parse_panel_entry(entry: str) -> Tuple[Optional[str], Optional[str]]:
    """'5b' -> ('5', 'b'), 'S6a' -> ('S6', 'a'), 'a' -> (None, 'a')"""
    entry = str(entry).strip()
    m = _PANEL_ENTRY_RE.match(entry)
    if m:
        return m.group(1), m.group(2).lower()
    # 纯字母如 'a'
    m = re.match(r"^([a-z])$", entry, re.I)
    if m:
        return None, m.group(1).lower()
    return None, None


def _collect_pairs_from_cols(
    group_df: pd.DataFrame,
    figure_col: str,
    panel_col: str,
) -> List[Tuple[str, str]]:
    """解析逗号分隔的 figure_id 和 panel_dwell 列，配对出 (figure_id, panel_letter)。

    figure_id 列: "Figure 5, Figure S10"
    panel 列:     "5b, 5e, S10c"
    -> [("Figure 5", "b"), ("Figure 5", "e"), ("Figure S10", "c")]
    """
    pairs: List[Tuple[str, str]] = []
    seen = set()

    if figure_col not in group_df.columns or panel_col not in group_df.columns:
        return pairs

    for _, row in group_df.iterrows():
        fid_val = row.get(figure_col)
        panel_val = row.get(panel_col)
        if fid_val is None or pd.isna(fid_val) or panel_val is None or pd.isna(panel_val):
            continue

        # 收集所有 figure_id 并建立 编号→完整ID 映射
        all_fids: List[str] = []
        fig_num_to_id: Dict[str, str] = {}
        for fid in str(fid_val).split(","):
            fid = fid.strip()
            if not fid:
                continue
            all_fids.append(fid)
            num = _extract_figure_number(fid)
            if num:
                fig_num_to_id[num.lower()] = fid

        if not all_fids:
            continue

        default_fid = all_fids[0]

        for entry in str(panel_val).split(","):
            entry = entry.strip()
            if not entry:
                continue
            fig_prefix, panel_letter = _parse_panel_entry(entry)
            if panel_letter is None:
                continue

            # 匹配 figure_prefix → 完整 figure_id
            if fig_prefix and fig_prefix.lower() in fig_num_to_id:
                fid = fig_num_to_id[fig_prefix.lower()]
            elif fig_prefix:
                # 宽松匹配
                fid = None
                for num, full in fig_num_to_id.items():
                    if fig_prefix.lower() == num:
                        fid = full
                        break
                if fid is None:
                    fid = default_fid
            else:
                fid = default_fid

            key = (fid, panel_letter)
            if key not in seen:
                seen.add(key)
                pairs.append(key)

    return pairs


def collect_candidate_panels(
    group_df: pd.DataFrame,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """返回 (dwell_current_pairs, dwell_time_pairs)。"""
    dc_pairs = _collect_pairs_from_cols(
        group_df, "figure_id_dwell_current", "panel_dwell_current",
    )
    dt_pairs = _collect_pairs_from_cols(
        group_df, "figure_id_dwell_time", "panel_dwell_time",
    )
    return dc_pairs, dt_pairs


# =============================================================================
# 7. 逐组处理
# =============================================================================


def evaluate_best_panel(
    vlm: DoubaoArkVLM,
    pdf_files: List[Path],
    candidates: List[Tuple[str, str]],
    group_info: Dict[str, Any],
    target_metric: str = "dwell_current",
    panel_loc_cache: Optional[Dict[Tuple[str, int], Dict[str, Tuple[float, float, float, float]]]] = None,
    figure_page_cache: Optional[Dict[str, Tuple[Optional[Path], Optional[int]]]] = None,
) -> Optional[Dict[str, Any]]:
    """对每组的所有候选面板逐一评估，选出最优。

    改进点：
    1. 先按 PDF 页面聚合 candidates，同一页只调用一次 VLM 识别 panel 位置。
    2. panel 坐标可跨 current/time、跨分组缓存复用。
    3. VLM 返回的 metric_type 先归一化，再判断是否为目标指标。
    """
    panel_loc_cache = panel_loc_cache if panel_loc_cache is not None else {}
    figure_page_cache = figure_page_cache if figure_page_cache is not None else {}

    evaluated: List[Dict[str, Any]] = []
    resolved_candidates: List[Dict[str, Any]] = []
    needed_panels_by_page: Dict[Tuple[str, int], set[str]] = defaultdict(set)
    pdf_by_page_key: Dict[Tuple[str, int], Path] = {}

    # 先定位所有候选 figure 的页面，并按页面聚合 panel。
    for figure_id, panel_label in candidates:
        panel_label = str(panel_label).strip().lower()
        print(f"    准备评估 [{target_metric}]: {figure_id} panel {panel_label}")

        fig_cache_key = str(figure_id).strip().lower()
        if fig_cache_key in figure_page_cache:
            pdf_path, page_no = figure_page_cache[fig_cache_key]
        else:
            pdf_path, page_no = find_figure_in_pdfs(pdf_files, figure_id)
            figure_page_cache[fig_cache_key] = (pdf_path, page_no)

        if not pdf_path or not page_no:
            print(f"      未找到 {figure_id} 在 PDF 中的位置，跳过")
            continue

        page_key = (str(pdf_path.resolve()), int(page_no))
        pdf_by_page_key[page_key] = pdf_path
        needed_panels_by_page[page_key].add(panel_label)
        resolved_candidates.append({
            "figure_id": figure_id,
            "panel": panel_label,
            "pdf_path": pdf_path,
            "page_no": int(page_no),
            "page_key": page_key,
        })

    if not resolved_candidates:
        print(f"    无可定位候选面板")
        return None

    # 同一 PDF 页面批量识别 panels，并缓存结果，避免重复调用 VLM。
    for page_key, target_panels in needed_panels_by_page.items():
        pdf_path = pdf_by_page_key[page_key]
        page_no = page_key[1]
        cached = panel_loc_cache.setdefault(page_key, {})
        missing_panels = sorted(p for p in target_panels if p not in cached)

        if missing_panels:
            print(
                f"      定位到: {pdf_path.name} 第 {page_no} 页；"
                f"批量识别 panels: {', '.join(missing_panels)}"
            )
            new_locs = vlm_identify_panels(vlm, pdf_path, page_no, missing_panels)
            cached.update(new_locs)
        else:
            print(f"      复用 panel 坐标缓存: {pdf_path.name} 第 {page_no} 页")

    # 再逐个 panel 裁剪并读取数值。
    for item in resolved_candidates:
        figure_id = item["figure_id"]
        panel_label = item["panel"]
        pdf_path = item["pdf_path"]
        page_no = item["page_no"]
        page_key = item["page_key"]

        print(f"    评估 [{target_metric}]: {figure_id} panel {panel_label}")
        panel_rect = panel_loc_cache.get(page_key, {}).get(panel_label)
        if not panel_rect:
            print(f"      未识别到 panel {panel_label} 位置，跳过")
            continue

        result = evaluate_panel(
            vlm, pdf_path, page_no, panel_rect,
            panel_label, figure_id, group_info,
            target_metric=target_metric,
        )
        if result is None:
            print(f"      评估失败，跳过")
            continue

        score = _safe_int(result.get("score", 0), default=0)
        metric_type_raw = result.get("metric_type", "")
        metric_type = normalize_metric_type(
            metric_type_raw,
            result.get("x_axis_label", ""),
            result.get("y_axis_label", ""),
        )
        is_target = metric_type == target_metric

        dc_val = result.get("dwell_current_value")
        dt_val = result.get("dwell_time_value")

        print(
            f"      metric_type_raw={metric_type_raw}, normalized={metric_type}, "
            f"score={score}, dc={dc_val}, dt={dt_val}"
        )

        evaluated.append({
            "figure_id": figure_id,
            "panel": panel_label,
            "pdf_name": pdf_path.name,
            "page_no": page_no,
            "metric_type": metric_type,
            "metric_type_raw": metric_type_raw,
            "score": score,
            "is_target": is_target,
            "dwell_current_value": dc_val,
            "dwell_current_unit": result.get("dwell_current_unit"),
            "dwell_time_value": dt_val,
            "dwell_time_unit": result.get("dwell_time_unit"),
            "chart_type": result.get("chart_type"),
            "x_axis_label": result.get("x_axis_label"),
            "y_axis_label": result.get("y_axis_label"),
            "selected_series_or_bar": result.get("selected_series_or_bar"),
            "value_note": result.get("value_note"),
            "reasoning": result.get("reasoning"),
        })

    if not evaluated:
        print(f"    无有效评估结果")
        return None

    # 只保留目标类型面板。
    target_panels = [e for e in evaluated if e["is_target"]]
    if not target_panels:
        print(f"    所有候选面板均非 {target_metric} 面板，无结果")
        return None

    # 优先选择已经读出原始值的 panel；若都没有值，则退回按 score 选最高。
    val_key = f"{target_metric}_value"
    panels_with_value = [e for e in target_panels if not _is_blank_value(e.get(val_key))]
    selection_pool = panels_with_value or target_panels
    best = max(selection_pool, key=lambda e: _safe_int(e.get("score", 0), default=0))
    print(f"    最优: {best['figure_id']} panel {best['panel']} (score={best['score']})")

    return {
        "best_figure_id": best["figure_id"],
        "best_panel": best["panel"],
        "best_pdf_name": best.get("pdf_name"),
        "best_page_no": best.get("page_no"),
        "best_score": best["score"],
        "best_value": best.get(val_key),
        "best_value_unit": best.get(f"{target_metric}_unit"),
        "best_chart_type": best.get("chart_type"),
        "candidates_evaluated": evaluated,
    }


# =============================================================================
# 7. 主处理
# =============================================================================

def process_paper_folder(
    folder: Path, vlm: DoubaoArkVLM,
) -> Optional[Path]:
    """处理单个 paper 文件夹。"""
    # 1. 找 filtered xlsx
    xlsx_files = list(folder.glob("*_filtered.xlsx"))
    if not xlsx_files:
        # 也尝试不带 filtered 后缀的
        xlsx_files = list(folder.glob("*_nanopore_blockade_signal_relations.xlsx"))
    if not xlsx_files:
        print(f"未找到 xlsx 文件于 {folder}")
        return None
    xlsx_path = xlsx_files[0]
    print(f"读取 Excel: {xlsx_path}")

    # 2. 找 PDF
    pdf_files = list(folder.glob("*.pdf"))
    if not pdf_files:
        print(f"未找到 PDF 文件于 {folder}")
        return None
    print(f"PDF 文件: {[p.name for p in pdf_files]}")

    # 3. 读取数据
    df = pd.read_excel(xlsx_path, sheet_name="signal_relation_records")
    print(f"总行数: {len(df)}")

    if df.empty:
        return None

    # 4. 按组处理
    existing_keys = [c for c in GROUP_KEY_COLS if c in df.columns]
    if not existing_keys:
        print("缺少分组列")
        return None

    groups = df.groupby(existing_keys, dropna=False)
    total_groups = len(groups)
    print(f"唯一分组数: {total_groups}")

    result_rows = []

    # 保留原始 df 的列结构，新增评估列
    output_df = df.copy()
    new_cols = [
        "best_fid_dwell_current", "best_panel_dwell_current",
        "best_dc_score", "best_dc_pdf", "best_dc_page",
        "vlm_dwell_current", "vlm_dwell_current_unit", "best_dc_chart_type",
        "best_fid_dwell_time", "best_panel_dwell_time",
        "best_dt_score", "best_dt_pdf", "best_dt_page",
        "vlm_dwell_time", "vlm_dwell_time_unit", "best_dt_chart_type",
        "candidate_evaluation_summary",
    ]
    for col in new_cols:
        output_df[col] = None

    # 缓存跨分组复用：
    # - figure_page_cache：避免反复在 PDF 文本中搜索同一个 Figure。
    # - panel_loc_cache：避免同一 PDF 页面重复调用 VLM 识别 panel 坐标。
    figure_page_cache: Dict[str, Tuple[Optional[Path], Optional[int]]] = {}
    panel_loc_cache: Dict[Tuple[str, int], Dict[str, Tuple[float, float, float, float]]] = {}

    for idx, (keys, group_df) in enumerate(groups, start=1):
        group_info = dict(zip(existing_keys, keys if isinstance(keys, tuple) else (keys,)))
        group_info["_context_summary"] = build_group_context(group_df, group_info)
        analyte = group_info.get("analyte", "?")
        print(f"\n[{idx}/{total_groups}] {analyte}")

        dc_pairs, dt_pairs = collect_candidate_panels(group_df)

        if not dc_pairs and not dt_pairs:
            print(f"  无候选面板，跳过")
            continue

        group_indices = group_df.index
        all_evaluated = []

        # ---- dwell_current ----
        if dc_pairs:
            print(f"  dwell_current 候选面板: {len(dc_pairs)} 个")
            for fid, panel in dc_pairs:
                print(f"    - {fid} panel {panel}")

            dc_result = evaluate_best_panel(
                vlm, pdf_files, dc_pairs, group_info,
                target_metric="dwell_current",
                panel_loc_cache=panel_loc_cache,
                figure_page_cache=figure_page_cache,
            )
            if dc_result:
                output_df.loc[group_indices, "best_fid_dwell_current"] = dc_result["best_figure_id"]
                output_df.loc[group_indices, "best_panel_dwell_current"] = dc_result["best_panel"]
                output_df.loc[group_indices, "best_dc_score"] = dc_result["best_score"]
                output_df.loc[group_indices, "best_dc_pdf"] = dc_result.get("best_pdf_name")
                output_df.loc[group_indices, "best_dc_page"] = dc_result.get("best_page_no")
                output_df.loc[group_indices, "vlm_dwell_current"] = dc_result["best_value"]
                output_df.loc[group_indices, "vlm_dwell_current_unit"] = dc_result.get("best_value_unit")
                output_df.loc[group_indices, "best_dc_chart_type"] = dc_result.get("best_chart_type")
                all_evaluated.extend(dc_result["candidates_evaluated"])

        # ---- dwell_time ----
        if dt_pairs:
            print(f"  dwell_time 候选面板: {len(dt_pairs)} 个")
            for fid, panel in dt_pairs:
                print(f"    - {fid} panel {panel}")

            dt_result = evaluate_best_panel(
                vlm, pdf_files, dt_pairs, group_info,
                target_metric="dwell_time",
                panel_loc_cache=panel_loc_cache,
                figure_page_cache=figure_page_cache,
            )
            if dt_result:
                output_df.loc[group_indices, "best_fid_dwell_time"] = dt_result["best_figure_id"]
                output_df.loc[group_indices, "best_panel_dwell_time"] = dt_result["best_panel"]
                output_df.loc[group_indices, "best_dt_score"] = dt_result["best_score"]
                output_df.loc[group_indices, "best_dt_pdf"] = dt_result.get("best_pdf_name")
                output_df.loc[group_indices, "best_dt_page"] = dt_result.get("best_page_no")
                output_df.loc[group_indices, "vlm_dwell_time"] = dt_result["best_value"]
                output_df.loc[group_indices, "vlm_dwell_time_unit"] = dt_result.get("best_value_unit")
                output_df.loc[group_indices, "best_dt_chart_type"] = dt_result.get("best_chart_type")
                all_evaluated.extend(dt_result["candidates_evaluated"])

        if all_evaluated:
            summary = json.dumps(all_evaluated, ensure_ascii=False)
            output_df.loc[group_indices, "candidate_evaluation_summary"] = summary

    # 5. 保存
    stem = xlsx_path.stem
    out_path = folder / f"{stem}_vlm_evaluated.xlsx"
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        output_df.to_excel(writer, index=False, sheet_name="signal_relation_records")
    print(f"\n已保存: {out_path}")
    return out_path


# =============================================================================
# 8. 文件查找
# =============================================================================

def find_paper_folders(root_dir: Path) -> List[Path]:
    """查找根目录下包含 filtered xlsx 和 PDF 的子文件夹。"""
    results = []
    for folder in sorted(p for p in root_dir.iterdir() if p.is_dir()):
        name = folder.name.strip()
        if not name or name in SKIP_FOLDERS:
            continue
        if name.startswith("."):
            continue
        if (name.startswith("__") and name.endswith("__")):
            continue
        # 检查是否包含需要的文件
        has_xlsx = bool(list(folder.glob("*_filtered.xlsx")) or
                        list(folder.glob("*_nanopore_blockade_signal_relations.xlsx")))
        has_pdf = bool(list(folder.glob("*.pdf")))
        if has_xlsx and has_pdf:
            results.append(folder)
    return results


# =============================================================================
# 9. main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="评估候选面板，选最优图读取 dwell_current 数据"
    )
    parser.add_argument(
        "root_dir", nargs="?", default=None,
        help="根目录路径（包含 paper 子文件夹）或直接指定 paper 文件夹",
    )
    parser.add_argument(
        "--folder", default=None,
        help="指定处理某个子文件夹",
    )
    parser.add_argument(
        "--vlm-model", default=None, help="VLM 模型 ID",
    )
    parser.add_argument(
        "--base-url", default=None, help="API 基础 URL",
    )
    args = parser.parse_args()

    vlm = DoubaoArkVLM(model=args.vlm_model, base_url=args.base_url)

    if not args.root_dir:
        root_dir = Path(__file__).resolve().parent
        print(f"未提供目录，默认使用脚本所在目录: {root_dir}")
    else:
        root_dir = Path(args.root_dir)

    if not root_dir.exists() or not root_dir.is_dir():
        print(f"目录不存在: {root_dir}")
        sys.exit(1)

    # 判断 root_dir 本身是否是 paper 文件夹（直接含 xlsx + pdf）
    has_xlsx = bool(list(root_dir.glob("*_filtered.xlsx")) or
                    list(root_dir.glob("*_nanopore_blockade_signal_relations.xlsx")))
    has_pdf = bool(list(root_dir.glob("*.pdf")))
    if has_xlsx and has_pdf:
        folders = [root_dir]
    else:
        folders = find_paper_folders(root_dir)

    if args.folder:
        folders = [f for f in folders if f.name == args.folder]
        if not folders:
            print(f"未找到文件夹: {args.folder}")
            sys.exit(1)

    if not folders:
        print("未找到包含 filtered xlsx 和 PDF 的文件夹")
        sys.exit(1)

    print(f"将处理 {len(folders)} 个文件夹: {[f.name for f in folders]}")

    for folder in folders:
        print(f"\n{'=' * 80}")
        print(f"处理文件夹: {folder.name}")
        print(f"{'=' * 80}")
        try:
            process_paper_folder(folder, vlm)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(f"处理失败: {exc}")

    print("\n全部完成")


if __name__ == "__main__":
    main()
 