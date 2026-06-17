"""Shared constants and helpers for the combine-insights-projects migration scripts."""

import csv

import boto3
import psycopg2

# ── Project UUIDs ─────────────────────────────────────────────────────────────
SET_1_ID = "37a0acdc-0a23-4a3b-ac4b-c47908b9ce17"
SET_2_ID = "722c3fbc-60f4-40f5-9d69-fdb0f3ca5421"
NEW_PROJECT_ID = "bb32aafc-2b75-4630-981f-77b92b1c2e43"
SOURCE_IDS = [SET_1_ID, SET_2_ID]

# ── AWS / DB constants ────────────────────────────────────────────────────────
REGION = "eu-west-1"
BUCKET = "biomage-originals-staging-242905224710"
DB_NAME = "aurora_db"
DB_USER = "dev_role"
DB_PORT = 5432


def get_connection(env, sandbox_id, aws_profile):
    """Open an IAM-authenticated psycopg2 connection via a local RDS tunnel."""
    print(f"Generating IAM token for {env}-{sandbox_id}...")
    session = boto3.Session(profile_name=aws_profile, region_name=REGION)
    rds = session.client("rds")
    response = rds.describe_db_cluster_endpoints(
        DBClusterIdentifier=f"aurora-cluster-{env}-{sandbox_id}",
        Filters=[{"Name": "db-cluster-endpoint-type", "Values": ["writer"]}],
    )
    endpoint = response["DBClusterEndpoints"][0]["Endpoint"]
    token = rds.generate_db_auth_token(endpoint, DB_PORT, DB_USER, REGION)
    print("Token generated, connecting...")
    return psycopg2.connect(
        host="localhost",
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=token,
        sslmode="require",
    )


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows)} rows → {path}")


def load_csv_map(path, key_fields, value_field):
    """Load a CSV into a dict keyed by tuple(key_fields) -> value_field."""
    result = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = tuple(row[k] for k in key_fields)
            result[key] = row[value_field]
    return result


def load_manifest(path):
    """Return dict: final_sample_name -> (experiment_id, original_sample_name)."""
    manifest = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            manifest[row["final_sample_name"]] = (
                row["experiment_id"],
                row["original_sample_name"],
            )
    return manifest
