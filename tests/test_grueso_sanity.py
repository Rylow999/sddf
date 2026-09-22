#!/usr/bin/env python3
"""Sanity check: spectrum_grueso.csv is consistent and what the paper claims."""
import csv
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "codigo"))
import numpy as np
from sddf_core import compute_spectral_curvature


def main():
    rows = list(csv.reader(open(
        Path(__file__).parent.parent / "datos/espectros/spectrum_grueso.csv")))
    k = np.array([float(r[0]) for r in rows[1:]])
    E = np.array([float(r[1]) for r in rows[1:]])
    G = compute_spectral_curvature(k, E)
    expected = 33.392638
    err = abs(G - expected) / expected
    assert err < 1e-6, f"G diff: {G} vs {expected} (err {err:.2e})"
    print(f"OK: G(spectrum_grueso) = {G:.6f} (expected {expected})")


if __name__ == "__main__":
    main()
