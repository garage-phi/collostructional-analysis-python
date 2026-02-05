"""
Test cases: Verification of the core module's basic structure.
"""

import pytest
import numpy as np
import pandas as pd

try:
    from core.collostructional_analysis import (
        AssociationStatsKernel,
        CollostructionalAnalysisMain,
    )
except ImportError:
    # If the core module cannot be imported, fail the tests
    pytest.fail("Could not import core module")

# It works without importing the validator from validation module
# from validation.validator_for_core import CollostructionalValidator


# =============================================================================
# 1. Normal Scenarios (Synthetic Data)
# Note: Tests without assets
# =============================================================================
class TestNormalScenarios:
    """Normal scenarios: Ensure calculations complete with minimal synthetic data"""

    def test_type1_simple_synthetic(self):
        """Type 1: Simple Collexeme Analysis runs with minimal data"""
        # Data creation
        df = pd.DataFrame(
            {
                "WORD": ["verb_a", "verb_b", "verb_c"],
                "FREQ_CORPUS": [1000, 2000, 500],
                "FREQ_CONST": [10, 5, 0],
            }
        )

        # Execution
        result = CollostructionalAnalysisMain.run(
            df,
            analysis_type=1,
            word_col="WORD",
            freq_corpus_col="FREQ_CORPUS",
            freq_const_col="FREQ_CONST",
            total_corpus_size=10000,
        )

        # Basic structure check
        assert not result.empty
        assert len(result) == 3
        expected_cols = ["WORD", "LLR", "FYE", "PMI", "LOGODDSRATIO"]
        for col in expected_cols:
            assert col in result.columns

        # Check data types of calculation results
        assert pd.api.types.is_float_dtype(result["LLR"])

    def test_type2_distinctive_synthetic(self):
        """Type 2: Distinctive Collexeme Analysis runs with minimal data"""
        df = pd.DataFrame(
            {
                "WORD": ["give", "give", "take", "take"],
                "CONST": [
                    "ditransitive",
                    "prepositional",
                    "ditransitive",
                    "prepositional",
                ],
            }
        )

        result = CollostructionalAnalysisMain.run(
            df, analysis_type=2, word_col="WORD", construction_col="CONST"
        )

        assert not result.empty
        # Distinctive analysis is pivoted,
        # so the number of rows should be the number of words (2)
        assert len(result) == 2
        assert "LLR" in result.columns

    def test_type2_comprehensive_scenarios(self):
        """
        Type 2: Distinctive Analysis with various input scenarios
        including raw tokens, frequency tables, and multiple constructions.
        """

        # ----------------------------------------------------------------
        # Case A: Raw Tokens input (2 Constructions)
        # Test the route where pivot_table is automatically applied
        # ----------------------------------------------------------------
        df_raw = pd.DataFrame(
            {
                "WORD": ["give", "give", "give", "take", "take"],
                "CXN": [
                    "ditransitive",
                    "ditransitive",
                    "prepositional",
                    "ditransitive",
                    "prepositional",
                ],
            }
        )
        # give: ditrans=2, prep=1
        # take: ditrans=1, prep=1

        res_raw = CollostructionalAnalysisMain.run(
            df_raw, analysis_type=2, word_col="WORD", construction_col="CXN"
        )

        assert not res_raw.empty
        assert len(res_raw) == 2  # words: give, take
        # Check if major metrics for 2 constructions are included
        assert "LOGODDSRATIO" in res_raw.columns
        assert "LLR" in res_raw.columns
        # Check if pivoted construction names are in columns
        assert "ditransitive" in res_raw.columns
        assert "prepositional" in res_raw.columns

        # ----------------------------------------------------------------
        # Case B: Frequency Table input (2 Constructions)
        # Test the route where user provides frequency columns directly
        # ----------------------------------------------------------------
        df_freq = pd.DataFrame(
            {"WORD": ["give", "take"], "FREQ_DITRANS": [100, 20], "FREQ_PREP": [50, 80]}
        )

        res_freq = CollostructionalAnalysisMain.run(
            df_freq,
            analysis_type=2,
            word_col="WORD",
            count_cols=["FREQ_DITRANS", "FREQ_PREP"],
        )

        assert not res_freq.empty
        # In non-raw mode, the provided column names are treated as construction names
        assert "FREQ_DITRANS" in res_freq.columns
        assert "Direction" in res_freq.columns

        # ----------------------------------------------------------------
        # Case C: Multiple Constructions (>2)
        # Test the route for handling multiple constructions
        # ----------------------------------------------------------------
        df_multi = pd.DataFrame(
            {
                "WORD": ["run", "walk"],
                "CXN_A": [10, 2],
                "CXN_B": [5, 10],
                "CXN_C": [2, 1],
            }
        )

        res_multi = CollostructionalAnalysisMain.run(
            df_multi,
            analysis_type=2,
            word_col="WORD",
            count_cols=["CXN_A", "CXN_B", "CXN_C"],
        )

        assert not res_multi.empty
        # Column Checks
        assert "SUMABSDEV" in res_multi.columns
        assert "LARGESTPREF" in res_multi.columns

        # LLR etc. are not output in multi-construction analysis (in this implementation),
        # so check that Residuals for each construction are included instead
        assert "CXN_A" in res_multi.columns  # Residual values are included

        # Check if LARGESTPREF is correctly selected
        # LargestPref should be CXN_A
        row_run = res_multi[res_multi["COLLOCATE"] == "run"].iloc[0]
        assert row_run["LARGESTPREF"] == "CXN_A"

    def test_type3_covarying_synthetic(self):
        """Type 3: Co-varying Collexeme Analysis runs with raw pair data"""
        # Data creation: List of co-occurring word pairs (Raw Tokens)
        # Cases:
        #  ("v1", "n1") -> 2 times
        #  ("v1", "n2")   -> 1 time
        #  ("v2", "n2")    -> 1 time
        #  ("v3", "n3")    -> 1 time
        df = pd.DataFrame(
            {
                "V_SLOT": ["v1", "v1", "v1", "v2", "v3"],
                "N_SLOT": ["n1", "n1", "n2", "n2", "n3"],
            }
        )

        # Total Corpus Size
        N = 100

        # Execution
        result = CollostructionalAnalysisMain.run(
            df,
            analysis_type=3,
            slot1_col="V_SLOT",
            slot2_col="N_SLOT",
            total_corpus_size=N,
        )

        # Check that the result is not empty
        assert not result.empty

        # Column Check
        # Expect columns after renaming as per rename_map in implementation
        expected_cols = ["V_SLOT", "N_SLOT", "Freq", "LLR", "PMI", "RELATION"]
        for col in expected_cols:
            assert col in result.columns

        # Verify aggregation logic
        # Find the pair "v1" - "n1"
        target = result[(result["V_SLOT"] == "v1") & (result["N_SLOT"] == "n1")].iloc[0]

        # Verify that Freq is 2 for ("v1", "n1")
        assert target["Freq"] == 2

        # Verify FREQOFSLOT1 value
        # FREQOFSLOT1: (a+b) = 3
        assert target["FREQOFSLOT1"] == 3
        # b = 3-a = 3-2 =1

        # Verify FREQOFSLOT2 value
        # FREQOFSLOT2: (a+c) = 2
        assert target["FREQOFSLOT2"] == 2
        # c = 2-a = 2-2 = 0

        # Check data types of calculation results
        res_signed = CollostructionalAnalysisMain.run(
            df,
            analysis_type=3,
            slot1_col="V_SLOT",
            slot2_col="N_SLOT",
            total_corpus_size=N,
            signed_metrics=True,
        )
        # Check that the "LLR" column is numeric type
        assert pd.api.types.is_numeric_dtype(res_signed["LLR"])

    def test_direct_calculation(self):
        """Direct calculation of association metrics from 2x2 table values"""
        # a=10, b=20, c=30, d=40
        res = AssociationStatsKernel.calculate_all_metrics(10, 20, 30, 40, 100)

        assert isinstance(res, pd.Series)
        assert "LLR" in res
        assert "PMI" in res
        # Check that LLR is not infinite or NaN
        assert np.isfinite(res["LLR"])
