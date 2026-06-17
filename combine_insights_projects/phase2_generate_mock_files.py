"""
Phase 2 — Generate mock sample files for upload to Insights UI.

Reads samples_and_files.csv (from Phase 1) to get sample names,
identifies duplicates across the two projects, applies _Set_1 / _Set_2
suffixes to duplicates only, and creates a folder per sample containing
minimal valid mock files:
    count_matrix.mtx  (matrixParse)
    all_genes.csv     (featuresParse)
    cell_metadata.csv (barcodesParse)

Output:
    mock_samples/<sample_name>/count_matrix.mtx
    mock_samples/<sample_name>/all_genes.csv
    mock_samples/<sample_name>/cell_metadata.csv
    sample_manifest.csv   — final sample names + source experiment
"""

import csv
import os
from collections import defaultdict

from utils import SET_1_ID, SET_2_ID

# Minimal valid Parse mock file contents
MOCK_MTX = """\
%%MatrixMarket matrix coordinate integer general
%
1 1 1
1 1 1
"""

MOCK_GENES = "MOCK_GENE,MOCK_GENE,Gene Expression\n"

MOCK_CELLS = "AAACCTGA,mock_sample\n"


def load_samples(csv_path):
    """Return a dict: experiment_id -> set of sample names (deduplicated per experiment)."""
    samples = defaultdict(set)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            samples[row["experiment_id"]].add(row["sample_name"])
    return samples


def build_final_names(samples_by_experiment):
    """
    Returns a list of (final_name, experiment_id, original_name) tuples.
    Duplicate names across experiments get _Set_1 / _Set_2 suffixes.
    """
    set1_names = samples_by_experiment.get(SET_1_ID, set())
    set2_names = samples_by_experiment.get(SET_2_ID, set())
    duplicates = set1_names & set2_names

    result = []
    for name in sorted(set1_names):
        final = f"{name}_Set_1" if name in duplicates else name
        result.append((final, SET_1_ID, name))
    for name in sorted(set2_names):
        final = f"{name}_Set_2" if name in duplicates else name
        result.append((final, SET_2_ID, name))
    return result


def create_mock_sample(out_dir, final_name):
    sample_dir = os.path.join(out_dir, final_name)
    os.makedirs(sample_dir, exist_ok=True)
    with open(os.path.join(sample_dir, "count_matrix.mtx"), "w") as f:
        f.write(MOCK_MTX)
    with open(os.path.join(sample_dir, "all_genes.csv"), "w") as f:
        f.write(MOCK_GENES)
    with open(os.path.join(sample_dir, "cell_metadata.csv"), "w") as f:
        f.write(MOCK_CELLS)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    samples_csv = os.path.join(script_dir, "samples_and_files.csv")
    out_dir = os.path.join(script_dir, "mock_samples")
    manifest_path = os.path.join(script_dir, "sample_manifest.csv")

    samples_by_experiment = load_samples(samples_csv)
    final_names = build_final_names(samples_by_experiment)

    print(f"Total samples: {len(final_names)}")
    duplicates = [n for n, _, orig in final_names if n.endswith("_Set_1") or n.endswith("_Set_2")]
    print(f"Renamed (duplicates): {len(duplicates)}")

    for final_name, exp_id, orig_name in final_names:
        create_mock_sample(out_dir, final_name)

    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["final_sample_name", "experiment_id", "original_sample_name"])
        writer.writeheader()
        for final_name, exp_id, orig_name in final_names:
            writer.writerow({
                "final_sample_name": final_name,
                "experiment_id": exp_id,
                "original_sample_name": orig_name,
            })

    print(f"Mock files written → {out_dir}/")
    print(f"Manifest written   → {manifest_path}")


if __name__ == "__main__":
    main()
