"""
Test Cases and Usage Examples for the Collostructional Analysis Methods

This file serves dual purposes:
1. Validation tests against Gries's v4.1 results
2. Usage examples for different analysis types:
    See the test cases below for practical usage patterns.

Requirements for Validation:
* Sample Data: Please download the sample input/output files 
    (e.g., 1.csv, 1_out.csv, 2a.csv, 2a_out.csv, etc.) from the official website of 
    Stefan Th. Gries: https://www.stgries.info/teaching/groningen/
* File Placement: Place these files in your current working directory 
    (typically the root of this repository where you run the script).
    * Note: The script utilizes `pd.read_csv("1.csv")`.
    * If you prefer to store data in a different directory, 
        please modify the file paths in the script accordingly.
* Environment: This script was primarily developed and verified on Google Colab.


Important Implementation Notes:

1. Log Odds Ratio Calculation:
   This script directly computes Log Odds Ratio based on the 2×2 contingency
   table definition: log((ad)/(bc)). Under this definition, Log Odds Ratio
   theoretically diverges in cases of perfect separation (b=0 or c=0).
   Prior research R implementations use glm(family = binomial), which may
   stop at finite values due to IRLS numerical convergence limits even
   when perfect separation occurs.

2. LLR (Log-Likelihood Ratio) Sign Convention:
   By default, this script returns absolute values for
   Log-Likelihood Ratio (LLR) and Fisher-Yates Exact test strength (FYE).
   Update (v1.1): Signed Metrics Mode You can optionally enable signed output by
   setting signed_metrics=True when calling the run() method.
   If enabled, both LLR and FYE will return negative values for repulsion patterns.
"""

import sys
import os
import dataclasses
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Union, Optional, Any

def _ensure_project_root_on_path():
    if "__file__" not in globals():
        # Notebook / Colab ?
        return

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

_ensure_project_root_on_path()

try:
    # 1. repository structure with package
    from core.collostructional_analysis import CollostructionalAnalysisMain
except ImportError:
    try:
        # 2. Loading in Colab or flat structure (when files are in the same location)
        from collostructional_analysis import CollostructionalAnalysisMain
    except ImportError as e:
        CollostructionalAnalysisMain = None
        _IMPORT_ERROR = e

def _require_analysis_main():
    if CollostructionalAnalysisMain is None:
        raise ImportError(
            "CollostructionalAnalysisMain could not be imported.\n"
            "Expected either:\n"
            "  - core/collostructional_analysis.py (package layout), or\n"
            "  - collostructional_analysis.py in the same directory.\n"
            f"Original error: {_IMPORT_ERROR}"
        )
    pass

# ==========================================
# 0. Path Handling & Imports
# ==========================================
# スタンドアロン実行時にプロジェクトルートが見えるようにする
if __name__ == "__main__":
    # このファイルが validation/ にある前提
    project_root = Path(__file__).parent
    if str(project_root) not in sys.path:
        sys.path.append(str(project_root))


# ==========================================
# 1. Dataclass / Result Structure
# ==========================================
@dataclasses.dataclass
class ValidationResult:
    """検証結果を保持するデータクラス"""
    success: bool
    errors: List[str]
    checked_rows: int
    
    def get_summary_message(self) -> str:
        if self.success:
            return "SUCCESS: All checks passed."
        else:
            return f"FAILURE: Found {len(self.errors)} types of mismatches."

# ==========================================
# 2. Core Validator Class
# ==========================================
class CollostructionalValidator:
    """
    Griesの実装結果とPython実装結果を比較する検証クラス。
    UI（print/pytest）から独立した純粋なロジックを持つ。
    """
    # Special Handling Column Names
    LLR_COLUMNS = {
        "LLR",
    }
    LOGODDS_COLUMNS = {
        "LOGODDSRATIO",
    }

    def __init__(self, atol: float = 1e-4, rtol: float = 1e-3):
        self.atol = atol
        self.rtol = rtol

    def validate(
        self, 
        result_df: pd.DataFrame, 
        expected_df: pd.DataFrame,
        key_col: Union[str, List[str]],
        col_mapping: Optional[Dict[str, str]] = None,
        exclude_cols: Optional[List[str]] = None,
        join_on: Optional[List[str]] = None
    ) -> ValidationResult:
        
        errors = []
        
        # --- Preprocessing ---
        res = result_df.copy()
        if col_mapping:
            res = res.rename(columns=col_mapping)
        expected = expected_df.copy()

        # print(f"len(res): {len(res)}, len(expected): {len(expected)}")
        # print(f"res.columns: {res.columns.tolist()}")
        # print(f"expected.columns: {expected.columns.tolist()}")

        # Join if requested (for Co-varying etc.)
        if join_on:
            expected = pd.merge(
                res[join_on].drop_duplicates(), 
                expected, 
                on=join_on, 
                how='inner'
            )
        # print(f"After join, len(expected): {len(expected)}")

        # Sort & Align
        #sort_keys = key_col if isinstance(key_col, list) else [key_col]
        # キーにNaNがあるとソートが不定になるため、最後尾に回す
        #res = res.sort_values(sort_keys, na_position="last")
        #expected = expected.sort_values(sort_keys, na_position="last")
        
        # print(f"indexes res: {res.index.tolist()}")
        # print(f"indexes expected: {expected.index.tolist()}")

        # 【重要修正】キーをインデックスに設定してアラインメントを保証する
        try:
            # 比較のために両方のDFのインデックスをキーに設定し、ソートする
            # これにより「行番号」ではなく「単語」で突き合わせが行われる
            res = res.set_index(key_col).sort_index()
            expected = expected.set_index(key_col).sort_index()
        except KeyError as e:
            return ValidationResult(False, [f"Key column missing or invalid: {e}"], 0)
        except Exception as e:
            return ValidationResult(False, [f"Index alignment failed: {e}"], 0)

        # Row Count Check
        if len(res) != len(expected):
            return ValidationResult(
                success=False, 
                errors=[f"Row count mismatch: Res={len(res)} vs Exp={len(expected)}"],
                checked_rows=0
            )

        # キーの集合が一致しているかチェック
        if not res.index.equals(expected.index):
             # 差分を表示
             missing_in_res = expected.index.difference(res.index)
             missing_in_exp = res.index.difference(expected.index)
             errs = []
             if not missing_in_res.empty: 
                 errs.append(f"Missing keys in Result (First 3): {missing_in_res[:3].tolist()}")
             if not missing_in_exp.empty: 
                 errs.append(f"Extra keys in Result (First 3): {missing_in_exp[:3].tolist()}")
             
             if errs:
                return ValidationResult(False, errs, len(res))
        
        # Column Intersection
        if exclude_cols:
            res = res.drop(columns=exclude_cols, errors='ignore')
            expected = expected.drop(columns=exclude_cols, errors='ignore')
        
        common_cols = [c for c in expected.columns if c in res.columns]
        if not common_cols:
            return ValidationResult(
                success=False, 
                errors=["No common columns found to compare."],
                checked_rows=len(res)
            )

        # --- Value Comparison ---
        for col in common_cols:
            r_vals = res[col]
            e_vals = expected[col]

            # Numeric Check
            if pd.api.types.is_numeric_dtype(r_vals) and pd.api.types.is_numeric_dtype(e_vals):
                msg = self._compare_numeric(r_vals, e_vals, col, res, key_col)
                if msg: errors.append(msg)
            
            # String Check
            else:
                msg = self._compare_string(r_vals, e_vals, col, res, key_col)
                if msg: errors.append(msg)

        return ValidationResult(
            success=(len(errors) == 0),
            errors=errors,
            checked_rows=len(res)
        )

    def _compare_numeric(self, r_vals, e_vals, col, df, key_col) -> Optional[str]:
        r_vals = r_vals.astype(float)
        e_vals = e_vals.astype(float)
        
        # Domain Logic: Signed Metrics (LLR) -> Compare Absolute
        if col in self.LLR_COLUMNS:
            is_close = np.isclose(np.abs(r_vals), np.abs(e_vals), rtol=self.rtol, atol=self.atol, equal_nan=True)
        
        # Domain Logic: LogOdds (Inf vs Finite Large)
        elif col in self.LOGODDS_COLUMNS:
            is_close_std = np.isclose(r_vals, e_vals, rtol=self.rtol, atol=self.atol, equal_nan=True)
            is_inf_match = np.isinf(r_vals) & (np.abs(e_vals) > 20) # Heuristic for R's glm limits
            is_close = is_close_std | is_inf_match
        
        else:
            is_close = np.isclose(r_vals, e_vals, rtol=self.rtol, atol=self.atol, equal_nan=True)

        if not is_close.all():
            diff_idx = np.where(~is_close)[0]
            details = self._format_error_details(diff_idx, r_vals, e_vals, df, key_col)
            return f"Numeric Mismatch in '{col}':\n{details}"
        return None

    def _compare_string(self, r_vals, e_vals, col, df, key_col) -> Optional[str]:
        try:
            r_str = r_vals.astype(str).str.lower().str.strip()
            e_str = e_vals.astype(str).str.lower().str.strip()
            is_eq = (r_str == e_str)
            if not is_eq.all():
                diff_idx = np.where(~is_eq)[0]
                details = self._format_error_details(diff_idx, r_vals, e_vals, df, key_col)
                return f"String Mismatch in '{col}':\n{details}"
        except Exception as e:
            if not r_vals.equals(e_vals):
                return f"Object Mismatch in '{col}': {e}"
        return None

    def _format_error_details(self, diff_idx, r_vals, e_vals, df, key_col) -> str:
        lines = []
        for i in diff_idx[:3]:
            # r_vals は Index が設定済みの Series なので .index[i] でキーが取れる
            idx_val = r_vals.index[i]
            val_res = r_vals.iloc[i]
            val_exp = e_vals.iloc[i]
            
            lines.append(f"    [Key={idx_val}] Res={val_res} != Exp={val_exp}")
            
        if len(diff_idx) > 3:
            lines.append(f"    ... and {len(diff_idx) - 3} more.")

        lines.append(f"total length: {len(df)}, key_col: {key_col}")
        return "\n".join(lines)


# ==========================================
# 3. CLI Runner (for Manual / Standalone)
# ==========================================
def run_cli():

    _require_analysis_main()

    print("=== Collostructional Analysis Validator CLI ===")

    # データディレクトリの自動推定 (repo/assets または repo/docs/assets)
    base_dir = Path(__file__).parent.parent
    assets_dir = base_dir / "assets"
    if not assets_dir.exists():
        assets_dir = base_dir / "assets" # fallback logic if needed
        # 今回の構成では assets が正

    print(f"Loading data from: {assets_dir}")
    if not assets_dir.exists():
        print("[Error] Assets directory not found.")
        return

    validator = CollostructionalValidator()

    # Define Cases to Run (Simplified for CLI)
    cases = [
        # --- Case 1: Simple Collexeme Analysis ---
        {
            "name": "Case 1 (Simple)",
            "file": "1.csv",
            "exp": "1_out.csv",
            "run_args": {
                "analysis_type": 1,
                "word_col": "WORD",
                "freq_corpus_col": "FREQ_WORD_in_CORPUS",
                "freq_const_col": "FREQ_WORD_in_DITRANSITIVE",
                "total_corpus_size": 138664,
            },
            "val_args": {
                "key_col": "WORD",
                "col_mapping": {
                    "a": "ditransitive",
                    "c": "OTHER",
                    "Direction": "RELATION",
                },
            },
        },
        # --- Case 2a: Distinctive (Raw Frequencies) ---
        {
            "name": "Case 2a (Distinctive Raw)",
            "file": "2a.csv",
            "exp": "2a_out.csv",
            "run_args": {
                "analysis_type": 2,
                "word_col": "Verb",
                "construction_col": "Construction",
            },
            "val_args": {
                "key_col": "WORD",
                "col_mapping": {
                    "Direction": "PREFERENCE",
                    "Verb": "WORD",
                    "a": "DITRANSITIVE",
                    "c": "PREP_DATIVE",
                },
            },
        },
        # --- Case 2b: Distinctive (Pre-calculated Frequencies) ---
        {
            "name": "Case 2b (Distinctive Freq)",
            "file": "2b.csv",
            "exp": "2b_out.csv",
            "run_args": {"analysis_type": 2, "word_col": "VERB"},
            "val_args": {
                "key_col": "WORD",
                "col_mapping": {
                    "Direction": "PREFERENCE",
                    "PMI": "MI",
                    "a": "DITRANSITIVE",
                    "c": "PREP_DATIVE",
                    "VERB": "WORD",
                },
            },
        },
        # --- Case 2c: Multiple Distinctive ---
        {
            "name": "Case 2c (Multiple)",
            "file": "2c.csv",
            "exp": "2c_out.csv",
            "run_args": {
                "analysis_type": 2,
                "word_col": "WORD",
                "construction_col": "CONSTRUCTION",
            },
            "val_args": {
                "key_col": "COLLOCATE",
                # No complex column mapping needed for 2c typically
            },
        },
        # --- Case 3: Co-varying Collexeme ---
        {
            "name": "Case 3 (Co-varying)",
            "file": "3.csv",
            "exp": "3_out.csv",
            "run_args": {"analysis_type": 3},
            "val_args": {
                "key_col": ["WORD_SLOT1", "WORD_SLOT2"],
                "col_mapping": {"Direction": "RELATION"},
                "join_on": ["WORD_SLOT1", "WORD_SLOT2"],
            },
        },
    ]

    for case in cases:
        print(f"\nRunning {case['name']}...")
        f_in = assets_dir / case['file']
        f_exp = assets_dir / case['exp']

        if not f_in.exists() or not f_exp.exists():
            print(f"  [Skip] Files not found: {case['file']}")
            continue

        try:
            df_in = pd.read_csv(f_in, sep='\t')
            df_exp = pd.read_csv(f_exp, sep='\t')

            # Run Analysis
            res_df = CollostructionalAnalysisMain.run(df_in, **case['run_args'])

            # Validate
            result = validator.validate(res_df, df_exp, **case['val_args'])

            if result.success:
                print(f"  [PASS] {result.get_summary_message()}")
            else:
                print(f"  [FAIL] {result.get_summary_message()}")
                for err in result.errors:
                    print(err)

        except Exception as e:
            print(f"  [Error] Execution failed: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    run_cli()
