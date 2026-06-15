"""
Phase 4 — Copy source files to destination keys in S3.

Reads copy_manifest.csv and performs server-side S3 copies
(no data leaves S3). Runs in parallel for speed.

Usage:
    python3 phase4_s3_copy.py [--profile <aws_profile>] [--workers 20]
"""

import argparse
import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3

BUCKET      = "biomage-originals-production-242905224710"
REGION      = "eu-west-1"
MAX_WORKERS = 20


def copy_one(s3, src, dst):
    s3.copy_object(
        CopySource={"Bucket": BUCKET, "Key": src},
        Bucket=BUCKET,
        Key=dst,
    )
    return src, dst


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", default="default")
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    manifest_path = os.path.join(script_dir, "copy_manifest.csv")

    with open(manifest_path) as f:
        rows = list(csv.DictReader(f))

    print(f"Copying {len(rows)} files → s3://{BUCKET}/")

    session = boto3.Session(profile_name=args.profile, region_name=REGION)
    s3 = session.client("s3")

    done, failed = 0, []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(copy_one, s3, row["src_s3_path"], row["dst_s3_path"]): row
            for row in rows
        }
        for future in as_completed(futures):
            row = futures[future]
            try:
                future.result()
                done += 1
                print(f"  [{done}/{len(rows)}] {row['final_sample_name']} / {row['sample_file_type']}")
            except Exception as e:
                failed.append((row, str(e)))
                print(f"  FAILED: {row['final_sample_name']} / {row['sample_file_type']}: {e}")

    print(f"\nDone: {done}/{len(rows)} copied.")
    if failed:
        print(f"Failed ({len(failed)}):")
        for row, err in failed:
            print(f"  {row['final_sample_name']} / {row['sample_file_type']}: {err}")


if __name__ == "__main__":
    main()
