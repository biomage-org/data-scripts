# Combine Insights Projects

## Goal

Combine two existing Trailmaker Insights projects into a single new project, without re-uploading count matrix files from scratch.

**Source projects:**
- Set_1: `b1322365-b917-4eda-b1e9-2964607efb58`
- Set_2: `ae905b9b-c6c4-4ac8-9ecb-be81b3ffde80`

**Target project:** `SO12129_Miranda_Lab_Set_1_and_2_Combined`

All samples from both projects are included. Samples whose names overlap across the two projects are renamed with a `_Set_1` / `_Set_2` suffix to avoid duplicates. A new metadata track called `Set` is added for all samples to record their project of origin (`Set_1` or `Set_2`). All other existing metadata tracks from both source projects are carried over.

---

## Plan

### Phase 1 — Reconnaissance (`phase1_reconnaissance.py`)
Query the production DB to collect, for both source projects:
- All sample IDs and names
- All `sampleFileId` S3 keys (one per file type per sample: `matrixParse`, `featuresParse`, `barcodesParse`)
- All metadata tracks and per-sample values

Outputs: `samples_and_files.csv`, `metadata.csv`

### Phase 2 — Create the new project (manual, Insights UI) (`phase2_generate_mock_files.py`)
Create a new empty project named `SO12129_Miranda_Lab_Set_1_and_2_Combined`.
Upload a minimal mock file triplet (`count_matrix.mtx`, `all_genes.csv`, `cell_metadata.csv`) for each sample, using the correct (possibly suffixed) sample names.

> **Why mock files?** The UI is the only way to create samples and have the system assign the new S3 UUIDs that Phase 3 captures. Uploading the real count matrices through the UI would be very slow for this many samples. Instead, minimal format-valid mock files are uploaded to trigger UUID assignment; Phase 4 then overwrites those S3 paths with the real data via a server-side copy — instant and free since no data leaves S3.

**Do not start the pipeline yet.**

### Phase 3 — Capture new UUIDs (`phase3_capture_new_uuids.py`)
Query the DB for the new project to get the mapping `sample_name → new sampleFileId` per file type.
Cross-reference with Phase 1 output to build a complete S3 copy manifest: `source_s3_path → destination_s3_path`.

Output: `copy_manifest.csv`

### Phase 4 — S3 copy (`phase4_s3_copy.sh`)
For each entry in the copy manifest, run:
```
aws s3 cp s3://biomage-originals-production-{accountId}/{src} \
           s3://biomage-originals-production-{accountId}/{dst}
```
No data leaves S3 — fast and free.

### Phase 5 — Insert metadata (`phase5_insert_metadata.py`)
Using the Phase 1 metadata and the sample ID mapping from Phase 3:
- Insert all metadata tracks from both source projects into the new experiment
- Insert per-sample values, respecting the renamed samples
- Insert the new `Set` track with `Set_1` / `Set_2` for all samples

### Phase 6 — Trigger pipeline & verify (manual)
Start the pipeline for the new project in Insights and confirm samples and metadata load correctly.

---

## Notes
- Each sample has 3 files in S3 (`matrixParse`, `featuresParse`, `barcodesParse`), each with its own UUID key in `biomage-originals`.
- Metadata is stored in PostgreSQL tables `metadata_track` and `sample_in_metadata_track_map`.
- Only overlapping sample names get the `_Set_1` / `_Set_2` suffix; non-overlapping names are kept as-is.
