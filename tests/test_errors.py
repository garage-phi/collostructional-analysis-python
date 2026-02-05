"""
Abnormal scenario tests for Collostructional Analysis core module.
This includes tests for error handling and edge cases.
"""

# import sys
# from pathlib import Path
import pandas as pd
import pytest

"""
# --- Path setup: see tests/conftest.py for details ---
def _ensure_project_root_on_path():
    if "__file__" not in globals(): return
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

_ensure_project_root_on_path()
"""

try:
    from core.collostructional_analysis import CollostructionalAnalysisMain
except ImportError:
    # Modify according to package structure
    pytest.fail("Could not import core module")

# importing the validator from validation module
from validation.validator_for_core import CollostructionalValidator


# --- Abnormal scenario tests ---

class TestAbnormalScenarios:
    """Abnormal scenarios and error handling tests"""

    def test_invalid_analysis_type(self):
        """
        Test that specifying a non-existent analysis type (e.g., 999) 
        raises a ValueError
        """
        df_dummy = pd.DataFrame({"A": [1], "B": [2]})
        
        # Expected exception: ValueError
        # The match argument specifies a part of the error message (can be a regex)
        with pytest.raises(ValueError, match="Analysis Type"):
            CollostructionalAnalysisMain.run(df_dummy, analysis_type=999)

    def test_empty_dataframe_is_rejected(self):
        """Test behavior when an empty DataFrame is passed"""
        df_empty = pd.DataFrame()
        
        # Implementation-dependent, but usually expect ValueError or similar
        with pytest.raises(ValueError):
            CollostructionalAnalysisMain.run(df_empty, analysis_type=1)

    # Type 1: Simple Collexeme Analysis
    def test_type1_missing_column(self):
        """
        Test that missing required columns raise appropriate errors.
        For example, if freq_const_col is missing.
        """
        df = pd.DataFrame({
            "WORD": ["give"],
            "FREQ_CORPUS": [100]
            # "FREQ_CONST" is missing
        })
        
        # KeyError or ValueError (adjust according to implementation)
        with pytest.raises(ValueError, match="FREQ_CONST"):
            CollostructionalAnalysisMain.run(
                df, 
                analysis_type=1,
                word_col="WORD",
                freq_corpus_col="FREQ_CORPUS",
                freq_const_col="FREQ_CONST", # This column is missing in df
                total_corpus_size=1000
            )

    def test_logodds_inf_with_small_expected_fails(self):
        df_res = pd.DataFrame({
            "WORD": ["a"],
            "LOGODDSRATIO": [float("inf")],
        })
        df_exp = pd.DataFrame({
            "WORD": ["a"],
            "LOGODDSRATIO": [5.0],  # too small
        })

        validator = CollostructionalValidator()

        # Execute validation
        result = validator.validate(df_res, df_exp, key_col="WORD")

        message = (
            "Should not succeed due to large value discrepancy"
        )

        # Assertion (Check that it "fails" as expected)
        assert result.success is False, message
        
        # Check that the error message contains appropriate content
        assert len(result.errors) > 0
        assert "Numeric Mismatch" in result.errors[0]
        assert "LOGODDSRATIO" in result.errors[0]   

    def test_logodds_inf_with_large_expected_passes(self):
        """
        [Reference] Conversely, when Expected is sufficiently large (>20),
        confirm that "Pass" occurs even with Inf and finite values.
        """
        df_res = pd.DataFrame({
            "WORD": ["a"],
            "LOGODDSRATIO": [float("inf")],
        })

        df_exp = pd.DataFrame({
            "WORD": ["a"],
            "LOGODDSRATIO": [25.0],          
        })
        # Since LOGODDSRATIO of df_exp exceeds 20, 
        # it should be considered equivalent to Inf

        validator = CollostructionalValidator()
        result = validator.validate(df_res, df_exp, key_col="WORD")

        message = "Should pass as expected value is sufficiently large"

        assert result.success is True, message

    def test_validator_key_mismatch_fails(self):
        # mismatch in key column

        df_res = pd.DataFrame({
            "WORD": ["a", "b"],
            "LLR": [1.0, 2.0],
        })
        df_exp = pd.DataFrame({
            "WORD": ["a", "c"],  # key mismatch
            "LLR": [1.0, 2.0],
        })

        validator = CollostructionalValidator()
        result = validator.validate(df_res, df_exp, key_col="WORD")

        message = "Should fail due to key column mismatch"

        assert result.success is False, message

    # Type 2: Distinctive Collexeme Analysis
    def test_type2_missing_column(self):
        """
        If the specified construction_col does not exist in the DataFrame,
        a ValueError should be raised.
        """
        df = pd.DataFrame({
            "Verb": ["give", "take"],
            # Not having "ConstCol" column
        })
        
        with pytest.raises(ValueError, match="ConstCol"):
            CollostructionalAnalysisMain.run(
                df,
                analysis_type=2,
                word_col="Verb",
                construction_col="ConstCol" # missing column
            )
    
    def test_type2_freq_mode_missing_count_col(self):
        """Type 2 (Freq Mode): When specified count_cols do not exist"""
        df = pd.DataFrame({
            "WORD": ["give", "take"],
            "CXN_A": [10, 20]
            # "CXN_B" is missing
        })

        # According to Pandas behavior, 
        # selecting a non-existent column raises a KeyError
        # (This error occurs inside the Analyzer's df[...] 
        # rather than in CollostructionalAnalysisMain)
        with pytest.raises(KeyError, match="CXN_B"):
            CollostructionalAnalysisMain.run(
                df,
                analysis_type=2,
                word_col="WORD",
                count_cols=["CXN_A", "CXN_B"] # missing columns
            )

    def test_type2_freq_mode_non_numeric_values(self):
        """
        Type 2 (Freq Mode): 
        If frequency columns contain non-numeric (non-convertible) values
        """
        df = pd.DataFrame({
            "WORD": ["give", "take"],
            "CXN_A": [10, 20],
            "CXN_B": ["many", "few"] # non-numeric (non-convertible) strings
        })

        # We expect a ValueError when .astype(int) is called internally
        with pytest.raises(ValueError):
            CollostructionalAnalysisMain.run(
                df,
                analysis_type=2,
                word_col="WORD",
                count_cols=["CXN_A", "CXN_B"]
            )

    # Type 3: Co-varying Collexeme Analysis
    def test_type3_missing_column(self):
        """Type 3: If the specified Slot column does not exist"""
        df = pd.DataFrame({
            "SlotA": ["a", "b"],
            "SlotB": ["x", "y"]
        })

        # "SlotA" exists but "WrongCol" does not
        with pytest.raises(ValueError, match="Specified column"):
            CollostructionalAnalysisMain.run(
                df,
                analysis_type=3,
                slot1_col="SlotA",
                slot2_col="WrongCol"
            )
    