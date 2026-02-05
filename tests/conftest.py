"""
Common configuration for pytest:
- Adds the project root to sys.path
- Provides a fixture for the assets directory
"""

import sys
from pathlib import Path
# import pandas as pd
import pytest

# ---------------------------------------------------------
# 1. Path Setup
# ---------------------------------------------------------
# Determine the project root, which is two levels up from tests/conftest.py
# .parent (tests folder) -> .parent (project root folder)
root_dir = Path(__file__).resolve().parent.parent

# Add the root to sys.path
# This allows "import core" and "import validation" to work from any test file
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))


# ---------------------------------------------------------
# 2. Fixtures (Shared resources for tests)
# ---------------------------------------------------------
@pytest.fixture(scope="session")
def data_dir():
    """
    Shared across all tests: 
    Fixture providing the path to the assets directory.
    
    Usage:
    def test_xxx(data_dir):
        csv_path = data_dir / "1.csv"

    Strategy: 
        To skip tests if the assets directory is 
        missing or contains no CSV files.
        Reason: Our Remote Repositry won't store csv files.
    """
    assets_path = root_dir / "assets"
    
    # 1. Check if assets directory exists
    if not assets_path.exists():
        pytest.skip("Skipping: Assets directory not found.")

    # 2. Check if there are any CSV files in the assets directory
    # Even if there is only .gitkeep, this will be False
    if not any(assets_path.glob("*.csv")):
        pytest.skip("Skipping: No CSV files found in assets directory.")

    return assets_path
