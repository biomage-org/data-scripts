"""
Phase 4 — Copy source files to destination keys in S3.

Reads copy_manifest.csv (or a custom --manifest file) and performs
server-side S3 copies (no data leaves S3). Runs in parallel for speed.

If any copies fail, the failed rows are written to failed.csv so they
can be retried with: python3 phase4_s3_copy.py --manifest failed.csv

Usage:
    python3 phase4_s3_copy.py [--profile <aws_profile>] [--workers 20] [--manifest copy_manifest.csv]
"""

import argparse
import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3

from utils import BUCKET, REGION

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
    parser.add_argument("--profile",  default="default")
    parser.add_argument("--workers",  type=int, default=MAX_WORKERS)
    parser.add_argument("--manifest", default=None, help="Path to manifest CSV (default: copy_manifest.csv next to this script)")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    manifest_path = args.manifest or os.path.join(script_dir, "copy_manifest.csv")

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

        failed_path = os.path.join(script_dir, "failed.csv")
        with open(failed_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["final_sample_name", "sample_file_type", "src_s3_path", "dst_s3_path"])
            writer.writeheader()
            writer.writerows(row for row, _ in failed)
        print(f"Failed entries written to {failed_path} — retry with: python3 phase4_s3_copy.py --manifest failed.csv")

        sys.exit(1)


if __name__ == "__main__":
    main()
