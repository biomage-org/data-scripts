"""
Phase 1 Reconnaissance — collect samples, S3 file keys, and metadata
for two source projects before combining them into a third.

Prerequisites:
    1. Open a tunnel to RDS in a separate terminal:
           biomage rds tunnel -i production -p <aws_profile>
    2. Install dependencies:
           pip install psycopg2-binary boto3

Usage:
    python phase1_reconnaissance.py [--env production] [--sandbox default] [--profile <aws_profile>]

Outputs (written to the same directory as this script):
    samples_and_files.csv   — sample IDs, names, file types, S3 keys
    metadata.csv            — metadata track values per sample
"""

import argparse
import os

import psycopg2.extras

from utils import SOURCE_IDS, get_connection, write_csv

SAMPLES_QUERY = """
SELECT
    s.experiment_id,
    s.id                AS sample_id,
    s.name              AS sample_name,
    sf.sample_file_type,
    sf.s3_path
FROM sample s
JOIN sample_to_sample_file_map m ON m.sample_id = s.id
JOIN sample_file sf              ON sf.id        = m.sample_file_id
WHERE s.experiment_id = ANY(%s::uuid[])
ORDER BY s.experiment_id, s.name, sf.sample_file_type;
"""

METADATA_QUERY = """
SELECT
    mt.experiment_id,
    sim.sample_id,
    s.name              AS sample_name,
    mt.key              AS metadata_track,
    sim.value
FROM metadata_track mt
JOIN sample_in_metadata_track_map sim ON sim.metadata_track_id = mt.id
JOIN sample s                         ON s.id                  = sim.sample_id
WHERE mt.experiment_id = ANY(%s::uuid[])
ORDER BY mt.experiment_id, s.name, mt.key;
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env",     default="production")
    parser.add_argument("--sandbox", default="default")
    parser.add_argument("--profile", default="default", help="AWS profile name")
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(__file__))
    conn = get_connection(args.env, args.sandbox, args.profile)

    with conn, conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        print("Querying samples and file keys...")
        cur.execute(SAMPLES_QUERY, (SOURCE_IDS,))
        samples = cur.fetchall()
        write_csv(
            os.path.join(out_dir, "samples_and_files.csv"),
            samples,
            ["experiment_id", "sample_id", "sample_name", "sample_file_type", "s3_path"],
        )

        print("Querying metadata tracks...")
        cur.execute(METADATA_QUERY, (SOURCE_IDS,))
        metadata = cur.fetchall()
        write_csv(
            os.path.join(out_dir, "metadata.csv"),
            metadata,
            ["experiment_id", "sample_id", "sample_name", "metadata_track", "value"],
        )

    conn.close()
    print("Done.")


if __name__ == "__main__":
    main()
