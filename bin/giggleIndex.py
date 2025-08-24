#!/usr/bin/env python

import argparse
from datetime import datetime
from giggle import Giggle

def main():
    parser = argparse.ArgumentParser(description="Create a Giggle index from BED files.")
    parser.add_argument('-i', '--inBed', required=True, help='Path to the input BED files directory')
    parser.add_argument('-o', '--outIdx', required=True, help='Path to the index output directory')
    args = parser.parse_args()

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] Creating Giggle index at {args.outIdx} from BED files in {args.inBed}")
    # Giggle.create('/home/xfd783/data/09_ALL_PUBLIC/database/unibind/v2022/mm10/test1', '/home/xfd783/data/09_ALL_PUBLIC/database/unibind/v2022/mm10/test/*.bed.gz')
    Giggle.create(args.outIdx, f"{args.inBed}/*.bed.gz")

if __name__ == "__main__":
    main()

# index not working
# https://github.com/ryanlayer/giggle/issues/21
