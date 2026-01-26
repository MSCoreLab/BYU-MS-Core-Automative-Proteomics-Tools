
from pathlib import Path
filename = "report.pg_matrix_57351_E25_20_2_350960_800.tsv"
loaded_stem = Path(filename).stem
print(f"Loaded Stem: '{loaded_stem}'")

double_stem = Path(loaded_stem).stem
print(f"Double Stem: '{double_stem}'")
