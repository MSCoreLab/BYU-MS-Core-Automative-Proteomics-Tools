import re
from pathlib import Path

filenames = [
    "report.pg_matrix_57349_E100_20_2_350960_600.tsv",
    "report.pg_matrix_57351_E25_20_2_350960_800.tsv",
    "report.pg_matrix_57353_E100_20_2_350960_800.tsv",
    "report.pg_matrix_57391_E25_30_4_440960_800.tsv",
    "report.pg_matrix_57392_E100_30_4_440960_800.tsv"
]

def get_mix_suffix(filename):
    # The EXACT regex from app.py
    match = re.search(r'(?:E[-_]?(?:25|100)|Y[-_]?(?:150|75))[-_](.*)', Path(filename).stem, re.IGNORECASE)
    return match.group(1) if match else None

print("Testing Regex Matches:")
for f in filenames:
    suffix = get_mix_suffix(f)
    print(f"File: {f}")
    print(f"Suffix: '{suffix}'")

print("\nSimulating Pairing:")
strict_pairs_dict = {}
singlets = []

# Mock lists
e25_samples = [f for f in filenames if "_E25_" in f]
e100_samples = [f for f in filenames if "_E100_" in f]

for s in e25_samples:
    suffix = get_mix_suffix(s)
    if suffix:
        if suffix not in strict_pairs_dict: strict_pairs_dict[suffix] = {}
        strict_pairs_dict[suffix]['E25'] = s
    else:
        print(f"FAILED to extract suffix for {s}")

for s in e100_samples:
    suffix = get_mix_suffix(s)
    if suffix:
        if suffix not in strict_pairs_dict: strict_pairs_dict[suffix] = {}
        strict_pairs_dict[suffix]['E100'] = s
    else:
        print(f"FAILED to extract suffix for {s}")

strict_pairs = []
for suffix, pair in strict_pairs_dict.items():
    if 'E25' in pair and 'E100' in pair:
        strict_pairs.append((pair['E25'], pair['E100']))
    else:
        singlets.extend(pair.values())

print(f"\nResults:")
print(f"Strict Pairs Found: {len(strict_pairs)}")
print(f"Singlets Found: {len(singlets)}")
for p in strict_pairs:
    print(f"Pair: {p}")
