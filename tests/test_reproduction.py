"""
Normal reproduction tests for Collostructional Analysis.
This verifies that the core analysis produces expected outputs
for given input datasets.
"""

# import sys
# from pathlib import Path
import pandas as pd
import pytest

# csvs: in (root)/assets 

"""
# Path setup: see tests/conftest.py for details
def _ensure_project_root_on_path():
    if "__file__" not in globals():
        # Notebook / Colab ?
        return

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))


_ensure_project_root_on_path()
"""

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

# importing the validator from validation module
from validation.validator_for_core import CollostructionalValidator


def load_data(data_dir, filename):
    p = data_dir / filename
    if not p.exists():
        pytest.skip(f"Missing file: {filename}")
    return pd.read_csv(p, sep="\t")


# Define test parameters
TEST_CASES = [
    pytest.param(
        "1.csv",
        "1_out.csv",
        {
            "analysis_type": 1,
            "word_col": "WORD",
            "freq_corpus_col": "FREQ_WORD_in_CORPUS",
            "freq_const_col": "FREQ_WORD_in_DITRANSITIVE",
            "total_corpus_size": 138664,
        },
        {
            "key_col": "WORD",
            "col_mapping": {"a": "ditransitive", "c": "OTHER", "Direction": "RELATION"},
        },
        id="case1",
    ),
    pytest.param(
        "2a.csv",
        "2a_out.csv",
        {"analysis_type": 2, "word_col": "Verb", "construction_col": "Construction"},
        {
            "key_col": "WORD",
            "col_mapping": {
                "Direction": "PREFERENCE",
                "Verb": "WORD",
                "a": "DITRANSITIVE",
                "c": "PREP_DATIVE",
            },
        },
        id="case2a",
    ),
    pytest.param(
        "2b.csv",
        "2b_out.csv",
        {"analysis_type": 2, "word_col": "VERB"},
        {
            "key_col": "WORD",
            "col_mapping": {
                "Direction": "PREFERENCE",
                "PMI": "MI",
                "a": "DITRANSITIVE",
                "c": "PREP_DATIVE",
                "VERB": "WORD",
            },
        },
        id="case2b",
    ),
    pytest.param(
        "2c.csv",
        "2c_out.csv",
        {"analysis_type": 2, "word_col": "WORD", "construction_col": "CONSTRUCTION"},
        {"key_col": "COLLOCATE"},
        id="case2c",
    ),
    # Case 3: Co-varying Collexeme
    pytest.param(
        "3.csv",
        "3_out.csv",
        {"analysis_type": 3},
        {
            "key_col": ["WORD_SLOT1", "WORD_SLOT2"],
            "col_mapping": {"Direction": "RELATION"},
            "join_on": ["WORD_SLOT1", "WORD_SLOT2"],
        },
        id="case3",
    ),
    # Other cases can be added similarly
]


@pytest.mark.parametrize("infile, expfile, run_kwargs, val_kwargs", TEST_CASES)
def test_reproduction(data_dir, infile, expfile, run_kwargs, val_kwargs):
    # 1. Preparation: Load input and expected output data
    df_in = load_data(data_dir, infile)
    df_exp = load_data(data_dir, expfile)

    # 2. Execution: Run the analysis
    result_df = CollostructionalAnalysisMain.run(df_in, **run_kwargs)

    # 3. Validation: Use the Validator class
    validator = CollostructionalValidator(atol=1e-4, rtol=1e-3)
    res_obj = validator.validate(result_df, df_exp, **val_kwargs)

    # 4. Decision: If the validation did not succeed, fail the test with error messages
    if not res_obj.success:
        # Concatenate error messages and fail the test
        error_msg = "\n\n".join(res_obj.errors)
        pytest.fail(f"Validation Mismatch:\n{error_msg}")
