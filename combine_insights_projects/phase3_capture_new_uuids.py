"""
Phase 3 — Capture new sampleFileId UUIDs and build S3 copy manifest.

Queries the new (combined) project from the DB to get the
sample_name → sampleFileId mapping per file type, then cross-references
with Phase 1 output to produce a copy manifest:
    src_s3_path  →  dst_s3_path  (per file type, per sample)

Prerequisites:
    - Tunnel open: biomage rds tunnel -i production -p <aws_profile>
    - Phase 1 outputs present: samples_and_files.csv, sample_manifest.csv

Usage:
    python3 phase3_capture_new_uuids.py [--env production] [--sandbox default] [--profile <aws_profile>]

Output:
    copy_manifest.csv — src_s3_path, dst_s3_path, sample_file_type, final_sample_name
"""

import argparse
import csv
import os
from collections import defaultdict

import boto3
import psycopg2
import psycopg2.extras

NEW_PROJECT_ID = "3b55335c-6f2b-409a-b854-adea40067095"

DB_NAME = "aurora_db"
DB_USER = "dev_role"
DB_PORT = 5432
REGION  = "eu-west-1"

NEW_PROJECT_QUERY = """
SELECT
    s.name              AS sample_name,
    sf.sample_file_type,
    sf.s3_path
FROM sample s
JOIN sample_to_sample_file_map m ON m.sample_id = s.id
JOIN sample_file sf              ON sf.id        = m.sample_file_id
WHERE s.experiment_id = %s::uuid
ORDER BY s.name, sf.sample_file_type;
"""


def get_connection(env, sandbox_id, aws_profile):
    session = boto3.Session(profile_name=aws_profile, region_name=REGION)
    rds = session.client("rds")
    response = rds.describe_db_cluster_endpoints(
        DBClusterIdentifier=f"aurora-cluster-{env}-{sandbox_id}",
        Filters=[{"Name": "db-cluster-endpoint-type", "Values": ["writer"]}],
    )
    endpoint = response["DBClusterEndpoints"][0]["Endpoint"]
    token = rds.generate_db_auth_token(endpoint, DB_PORT, DB_USER, REGION)
    return psycopg2.connect(
        host="localhost", port=DB_PORT, dbname=DB_NAME,
        user=DB_USER, password=token, sslmode="require",
    )


def load_csv(path, key_fields, value_field):
    """Load a CSV into a dict keyed by tuple(key_fields) -> value_field."""
    result = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = tuple(row[k] for k in key_fields)
            result[key] = row[value_field]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env",     default="production")
    parser.add_argument("--sandbox", default="default")
    parser.add_argument("--profile", default="default")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Load Phase 1 outputs
    # samples_and_files: (experiment_id, original_sample_name, file_type) -> src_s3_path
    src_keys = load_csv(
        os.path.join(script_dir, "samples_and_files.csv"),
        key_fields=["experiment_id", "sample_name", "sample_file_type"],
        value_field="s3_path",
    )

    # sample_manifest: final_sample_name -> (experiment_id, original_sample_name)
    manifest = {}
    with open(os.path.join(script_dir, "sample_manifest.csv")) as f:
        for row in csv.DictReader(f):
            manifest[row["final_sample_name"]] = (row["experiment_id"], row["original_sample_name"])

    # Query new project
    print("Querying new project for destination UUIDs...")
    conn = get_connection(args.env, args.sandbox, args.profile)
    with conn, conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(NEW_PROJECT_QUERY, (NEW_PROJECT_ID,))
        new_rows = cur.fetchall()
    conn.close()
    print(f"  Found {len(new_rows)} file records in new project")

    # Build copy manifest
    manifest_rows = []
    missing = []

    for row in new_rows:
        final_name = row["sample_name"]
        file_type  = row["sample_file_type"]
        dst_path   = row["s3_path"]

        if final_name not in manifest:
            missing.append(f"  No manifest entry for sample: {final_name}")
            continue

        exp_id, orig_name = manifest[final_name]
        src_path = src_keys.get((exp_id, orig_name, file_type))

        if not src_path:
            missing.append(f"  No source file for: {final_name} / {file_type}")
            continue

        manifest_rows.append({
            "final_sample_name": final_name,
            "sample_file_type":  file_type,
            "src_s3_path":       src_path,
            "dst_s3_path":       dst_path,
        })

    if missing:
        print("WARNINGS — could not resolve the following:")
        for m in missing:
            print(m)

    out_path = os.path.join(script_dir, "copy_manifest.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["final_sample_name", "sample_file_type", "src_s3_path", "dst_s3_path"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Copy manifest written → {out_path}  ({len(manifest_rows)} entries)")


if __name__ == "__main__":
    main()
