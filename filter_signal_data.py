#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_signal_data.py

过滤已提取的纳米孔阻断信号数据：
1. 保留 analyte_concentration_value 为 5、10、15、20 的行
2. 保留 voltage_mV 为 20 mV 倍数的行

用法：
    python filter_signal_data.py "D:/Code/Python/extract_peptide"

    # 处理单个 Excel 文件
    python filter_signal_data.py --excel "D:/path/to/file.xlsx"
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

VALID_CONCENTRATIONS = {5, 10, 15, 20}
VOLTAGE_STEP = 20

_FIG_NUM_RE = re.compile(r"(?:Figure|Fig\.?)\s*(S?\d+)", re.I)
_PANEL_ENTRY_RE = re.compile(r"^(s?\d+)([a-z])$", re.I)


def should_skip_folder(folder: Path) -> bool:
    name = folder.name.strip()
    if not name:
        return True
    lower = name.lower()
    if lower in {"__pycache__", ".git", ".idea", ".vscode", ".ipynb_checkpoints"}:
        return True
    if name.startswith(".") or (name.startswith("__") and name.endswith("__")):
        return True
    return False


def find_excel_files(root_dir: Path) -> List[Tuple[Path, Path]]:
    """查找所有子文件夹中的 *_nanopore_blockade_signal_relations.xlsx。"""
    results = []
    for folder in sorted(p for p in root_dir.iterdir() if p.is_dir()):
        if should_skip_folder(folder):
            continue
        for excel_path in folder.glob("*_nanopore_blockade_signal_relations.xlsx"):
            results.append((folder, excel_path))
    return results


def extract_voltage_mv(row: pd.Series) -> Optional[float]:
    """尝试从行中提取电压数值 (mV)。"""
    vm = row.get("voltage_mV")
    if vm is not None and pd.notna(vm):
        try:
            return abs(float(vm))
        except (ValueError, TypeError):
            pass

    vr = row.get("voltage_raw")
    if vr is not None and pd.notna(vr):
        m = re.search(r"[+-]?(\d+(?:\.\d+)?)\s*mV", str(vr))
        if m:
            try:
                return abs(float(m.group(1)))
            except (ValueError, TypeError):
                pass
    return None


def _extract_figure_number(figure_id: str) -> Optional[str]:
    """'Figure 5' -> '5', 'Figure S10' -> 'S10'"""
    m = _FIG_NUM_RE.search(str(figure_id))
    return m.group(1) if m else None


def _extract_panel_fig_prefix(panel_entry: str) -> Optional[str]:
    """'5b' -> '5', 'S6a' -> 'S6', 'a' -> None"""
    m = _PANEL_ENTRY_RE.match(str(panel_entry).strip())
    return m.group(1) if m else None


def _clean_panel_column(panel_text: object) -> Optional[str]:
    """清理 panel 列的括号描述，只保留出处如 S6a。

    '5b (Ib/I0 histogram subpanel)' -> '5b'
    'S6a, 5b (tD histogram subpanel)' -> 'S6a, 5b'
    """
    if panel_text is None or pd.isna(panel_text):
        return None
    parts = []
    for entry in str(panel_text).split(","):
        entry = re.sub(r"\s*\([^)]*\)\s*", "", entry).strip()
        if entry:
            parts.append(entry)
    return ", ".join(parts) if parts else None


def _match_figures_to_panels(figure_text: object, panel_text: object) -> Optional[str]:
    """删掉 figure_id 列中在 panel 列里没有对应条目的 figure。

    figure_id_dwell_current: "Figure 5, Figure S10, Figure S11"
    panel_dwell_current:     "5b, 5e, 5f, 5g"
    -> "Figure 5"  （S10/S11 没有对应 panel，删掉）
    """
    if figure_text is None or pd.isna(figure_text):
        return None
    if panel_text is None or pd.isna(panel_text):
        return str(figure_text).strip() or None

    panel_fig_nums: set[str] = set()
    for entry in str(panel_text).split(","):
        prefix = _extract_panel_fig_prefix(entry.strip())
        if prefix:
            panel_fig_nums.add(prefix.lower())

    if not panel_fig_nums:
        return str(figure_text).strip() or None

    kept: list[str] = []
    for fid in str(figure_text).split(","):
        fid = fid.strip()
        if not fid:
            continue
        num = _extract_figure_number(fid)
        if num and num.lower() in panel_fig_nums:
            kept.append(fid)

    return ", ".join(kept) if kept else None


def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """过滤 DataFrame：浓度限于 5/10/15/20，电压为 20 mV 倍数。"""
    if df.empty:
        return df

    mask = pd.Series(True, index=df.index)

    # 1. 浓度过滤
    if "analyte_concentration_value" in df.columns:
        conc_mask = df["analyte_concentration_value"].apply(
            lambda x: x in VALID_CONCENTRATIONS if pd.notna(x) else False
        )
        mask = mask & conc_mask
    else:
        print("  警告: 未找到 analyte_concentration_value 列，跳过浓度过滤")

    # 2. 电压过滤：20 mV 的倍数
    voltage_mv = df.apply(extract_voltage_mv, axis=1)
    voltage_mask = voltage_mv.apply(
        lambda v: v is not None and abs(v) > 0 and abs(v) % VOLTAGE_STEP == 0
    )
    mask = mask & voltage_mask

    df = df[mask].copy()

    # ---- 清理 panel 列：去掉括号描述，只保留出处如 S6a ----
    for panel_col in ("panel_dwell_current", "panel_dwell_time"):
        if panel_col in df.columns:
            df[panel_col] = df[panel_col].apply(_clean_panel_column)

    # ---- 匹配 figure_id_dwell_* 与 panel_dwell_*，删掉没有对应 panel 的 figure ----
    for fig_col, panel_col in [
        ("figure_id_dwell_current", "panel_dwell_current"),
        ("figure_id_dwell_time", "panel_dwell_time"),
    ]:
        if fig_col in df.columns and panel_col in df.columns:
            df[fig_col] = df.apply(
                lambda row, fc=fig_col, pc=panel_col: _match_figures_to_panels(row.get(fc), row.get(pc)),
                axis=1,
            )

    return df


def process_excel(excel_path: Path, folder_path: Path) -> bool:
    print(f"\n处理: {excel_path}")

    df = pd.read_excel(excel_path, sheet_name="signal_relation_records")
    original_count = len(df)
    print(f"  原始行数: {original_count}")

    if df.empty:
        return False

    filtered = filter_dataframe(df)
    filtered_count = len(filtered)
    removed = original_count - filtered_count
    print(f"  过滤后行数: {filtered_count}（移除 {removed} 行）")

    output_path = folder_path / f"{excel_path.stem}_filtered.xlsx"
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        filtered.to_excel(writer, index=False, sheet_name="signal_relation_records")
    print(f"  已保存: {output_path}")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="过滤纳米孔阻断信号数据")
    parser.add_argument("root_dir", nargs="?", default=None, help="根目录路径")
    parser.add_argument("--excel", default=None, help="单个 Excel 文件路径")
    args = parser.parse_args()

    if args.excel:
        excel_path = Path(args.excel)
        if not excel_path.exists():
            print(f"文件不存在: {excel_path}")
            sys.exit(1)
        process_excel(excel_path, excel_path.parent)
        return

    if args.root_dir:
        root_dir = Path(args.root_dir)
    else:
        root_dir = Path(__file__).resolve().parent
        print(f"未提供根目录，默认使用脚本所在目录: {root_dir}")

    if not root_dir.exists() or not root_dir.is_dir():
        print(f"根目录不存在: {root_dir}")
        sys.exit(1)

    excel_files = find_excel_files(root_dir)
    if not excel_files:
        print(f"未找到任何 Excel 文件")
        sys.exit(1)

    print(f"找到 {len(excel_files)} 个 Excel 文件")
    success = 0
    for folder_path, excel_path in excel_files:
        if process_excel(excel_path, folder_path):
            success += 1

    print(f"\n处理完成: {success}/{len(excel_files)} 个文件")


if __name__ == "__main__":
    main()
