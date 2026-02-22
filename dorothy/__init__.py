"""
Dorothy - Crystallography Teaching Tool

Honoring Dorothy Hodgkin's pioneering work in X-ray crystallography.
"""

import sys
from pathlib import Path

__version__ = "0.6.0-beta"
__author__ = "Your Name"


def _base_dir() -> Path:
    """Project root — works both in dev and PyInstaller frozen mode."""
    if getattr(sys, 'frozen', False):
        return Path(sys._MEIPASS)
    return Path(__file__).parent.parent
