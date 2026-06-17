"""
Phase 5 — Insert sample-level metadata into the new combined project.

For each metadata track from the two source projects:
  - Creates a metadata_track row in the new experiment
  - Inserts per-sample values (N.A. for samples whose source project
    didn't have that track)

Also inserts a new "Set" track: Set_1 / Set_2 based on origin project.

Prerequisites:
    - Tunnel open: biomage rds tunnel -i production -p <aws_profile>
    - Phase 1 outputs: metadata.csv, sample_manifest.csv
    - Phase 3 output: copy_manifest.csv (to confirm sample count)

Usage:
    python3 phase5_insert_metadata.py [--env production] [--sandbox default] [--profile <aws_profile>] [--dry-run]

    Use --dry-run to print what would be inserted without touching the DB.
"""

import argparse
import os

import psycopg2.extras

from utils import NEW_PROJECT_ID, SET_1_ID, get_connection, load_csv_map, load_manifest

NEW_SAMPLES_QUERY = """
SELECT id AS sample_id, name AS sample_name
FROM sample
WHERE experiment_id = %s::uuid
ORDER BY name;
"""


def get_all_track_keys(metadata):
    """Return sorted list of unique metadata track keys across both projects."""
    return sorted({track for (_, _, track) in metadata.keys()})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env",     default="production")
    parser.add_argument("--sandbox", default="default")
    parser.add_argument("--profile", default="default")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without writing to DB")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    metadata   = load_csv_map(os.path.join(script_dir, "metadata.csv"), ["experiment_id", "sample_name", "metadata_track"], "value")
    manifest   = load_manifest(os.path.join(script_dir, "sample_manifest.csv"))
    track_keys = get_all_track_keys(metadata)

    # Add the new Set track
    all_tracks = track_keys + ["Set"]
    print(f"Tracks to insert: {all_tracks}")

    conn = get_connection(args.env, args.sandbox, args.profile)
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        # Fetch new project's samples
        cur.execute(NEW_SAMPLES_QUERY, (NEW_PROJECT_ID,))
        new_samples = {row["sample_name"]: row["sample_id"] for row in cur.fetchall()}

    print(f"Samples in new project: {len(new_samples)}")

    if args.dry_run:
        print("\n--- DRY RUN --- (no DB writes)\n")

    try:
        with conn:
            with conn.cursor() as cur:
                for track_key in all_tracks:
                    if args.dry_run:
                        print(f"INSERT metadata_track: experiment={NEW_PROJECT_ID}, key={track_key}")
                    else:
                        cur.execute(
                            "INSERT INTO metadata_track (experiment_id, key) VALUES (%s::uuid, %s) RETURNING id;",
                            (NEW_PROJECT_ID, track_key),
                        )
                    track_id = cur.fetchone()[0] if not args.dry_run else None

                    rows_inserted = 0
                    for final_name, new_sample_id in new_samples.items():
                        if track_key == "Set":
                            exp_id, _ = manifest.get(final_name, (None, None))
                            value = "Set_1" if exp_id == SET_1_ID else "Set_2"
                        else:
                            exp_id, orig_name = manifest.get(final_name, (None, None))
                            value = metadata.get((exp_id, orig_name, track_key), "N.A.")

                        if args.dry_run:
                            print(f"  INSERT map: sample={final_name}, track={track_key}, value={value}")
                        else:
                            cur.execute(
                                """INSERT INTO sample_in_metadata_track_map
                                       (metadata_track_id, sample_id, value)
                                   VALUES (%s, %s::uuid, %s);""",
                                (track_id, new_sample_id, value),
                            )
                        rows_inserted += 1

                    if not args.dry_run:
                        print(f"  Inserted track '{track_key}' with {rows_inserted} sample values")

    except Exception as e:
        print(f"ERROR: {e}")
        conn.rollback()
        raise
    finally:
        conn.close()

    print("\nDone." if not args.dry_run else "\nDry run complete — no changes made.")


if __name__ == "__main__":
    main()
