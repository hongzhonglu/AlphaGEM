"""Compatibility entry point; use the same pipeline as AlphaGEM.py."""

import sys
from pathlib import Path

_repository_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repository_root))
from AlphaGEM import main


if __name__ == "__main__":
    # Preserve the historical spellings while sharing defaults and solver setup.
    sys.argv[1:] = [
        "structure_alignment" if arg == "structure alignment" else
        "--esp" if arg == "--spe" else arg
        for arg in sys.argv[1:]
    ]
    main()
