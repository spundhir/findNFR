#!/usr/bin/env python

import argparse
import os
import sys
import tempfile
import subprocess
import gzip
from datetime import datetime

def usage():
    print("""Program: giggleAna (annotate genomic regions for overlap with giggle database)
Author: BRIC, University of Copenhagen, Denmark
Version: 1.0
Contact: pundhir@binf.ku.dk
Usage: giggleAna.py -i <file> [OPTIONS]
Options:
 -i <file>   [input file having genomic regions in BED format (or stdin)]
[OPTIONS]
 -g <string> [genome for which to perform the analysis (mm10, hg38, danRer7; default: mm10)]
 -d <string> [database to query (atac, jaspar, gtrd, unibind, homer; default: unibind)]
 -h          [help]
""")
    sys.exit(0)
    
def process_bed(input_bed, output_bed_gz, keep_id):
    counter = 1
    with open(input_bed, 'rt') as fin, gzip.open(output_bed_gz, 'wt') as fout:
        for line in fin:
            fields = line.rstrip('\n').split('\t')
            # Ensure at least 6 columns
            while len(fields) < 6:
                fields.append('.')
            if fields[5] not in ['+', '-']:
                fields[5] = '.'
            if keep_id:
                out_fields = [fields[0], fields[1], fields[2], fields[3], '1', fields[5]]
            else:
                out_fields = [fields[0], fields[1], fields[2], f'REGION{counter}', '1', fields[5]]
            fout.write('\t'.join(out_fields) + '\n')
            counter += 1
            
def append_total_regions(giggle_output, target_total, output_file):
    with open(giggle_output, 'rt') as fin, open(output_file, 'wt') as fout:
        for line in fin:
            line = line.rstrip('\n')
            if 'odds_ratio' in line:
                fout.write(f"{line}\ttotal_regions\n")
            else:
                fout.write(f"{line}\t{target_total}\n")
                
def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', dest='bedfile', required=False)
    parser.add_argument('-g', dest='genome', default='mm10')
    parser.add_argument('-d', dest='database', default='unibind')
    parser.add_argument('-h', dest='help', action='store_true')
    args = parser.parse_args()

    if args.help or not args.bedfile:
        usage()

    # Set directories based on genome
    genome = args.genome
    database = args.database

    print(f"Populating files based on input genome, {genome} ({datetime.now().strftime('%c')}).. ", file=sys.stderr, end='')
    dirs = {
        "mm10": {
            "atac": "/home/xfd783/data/09_ALL_PUBLIC/database/chipAtlas/atac/mm10/giggle/",
            "jaspar": "/home/xfd783/data/09_ALL_PUBLIC/database/jaspar/mm10/giggle/",
            "gtrd": "/home/xfd783/data/09_ALL_PUBLIC/database/gtrd/v21.12/mm10/giggle/",
            "unibind": "/home/xfd783/data/09_ALL_PUBLIC/database/unibind/v2022/mm10/giggle/",
            "homer": "/home/xfd783/data/09_ALL_PUBLIC/database/homer/mm10/giggle"
        },
        "hg38": {
            "atac": "/home/xfd783/data/09_ALL_PUBLIC/database/chipAtlas/atac/hg38/giggle/",
            "jaspar": "/home/xfd783/data/09_ALL_PUBLIC/database/jaspar/hg38/giggle/",
            "gtrd": "/home/xfd783/data/09_ALL_PUBLIC/database/gtrd/v21.12/hg38/giggle/",
            "unibind": "/home/xfd783/data/09_ALL_PUBLIC/database/unibind/v2022/hg38/giggle/",
            "homer": "/home/xfd783/data/09_ALL_PUBLIC/database/homer/hg38/giggle"
        },
        "danRer7": {
            "atac": "",
            "jaspar": "",
            "gtrd": "",
            "unibind": "",
            "homer": ""
        }
    }

    if genome not in dirs:
        print("Presently the program only supports analysis for mm10 or hg38\n")
        usage()

    db_dir = dirs[genome].get(database, "")
    if not db_dir:
        print("Invalid database or genome combination.")
        usage()
    print("done", file=sys.stderr)

    # Create region file depending on if the input is from file or STDIN
    print("Create region file depending on if the input is from file or STDIN... ", file=sys.stderr, end='')
    with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as tmp_bed:
        tmp_bed_name = tmp_bed.name
        if args.bedfile == "stdin":
            for line in sys.stdin:
                tmp_bed.write(line.encode())
        elif os.path.isfile(args.bedfile):
            if args.bedfile.endswith('.gz'):
                subprocess.run(['zless', args.bedfile], stdout=tmp_bed)
            else:
                with open(args.bedfile, 'rb') as f:
                    tmp_bed.write(f.read())
        else:
            usage()
    print("done", file=sys.stderr)

    # Count rows and unique IDs
    row_count = int(subprocess.check_output(f"zless {tmp_bed_name} | wc -l", shell=True).decode().strip())
    id_count = int(subprocess.check_output(f"zless {tmp_bed_name} | cut -f 4 | sort | uniq | wc -l", shell=True).decode().strip())
    keep_id = 1 if row_count == id_count else 0

    # Ensure peak ids are unique and compress
    tmp_bed_gz = tmp_bed_name + ".gz"
    perl_cmd = (
        f"zless {tmp_bed_name} | "
        "perl -ane 'BEGIN { $counter=1; } "
        "if($F[5]!~/^\\+$/ && $F[5]!~/^\\-$/) { $F[5]=\".\"; } "
        f"if({keep_id}) {{ print \"$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t1\\t$F[5]\"; }} "
        "else { print \"$F[0]\\t$F[1]\\t$F[2]\\tREGION$counter\\t1\\t$F[5]\"; } "
        "print \"\\n\"; $counter++;' | bgzip > {tmp_bed_gz}"
    )
    subprocess.run(perl_cmd, shell=True, executable='/bin/bash')

    # Perform giggle enrichment analysis
    target_total = int(subprocess.check_output(f"zcat {tmp_bed_gz} | wc -l", shell=True).decode().strip())
    giggle_cmd = ""
    if os.path.isdir(db_dir):
        giggle_cmd = f"giggle search -i {db_dir} -q {tmp_bed_gz} -s"
        if database == "atac":
            giggle_cmd += " -f 05"
        giggle_cmd += (
            " | perl -ane 'chomp($_); "
            "if($_=~/odds_ratio/) { print \"$_\\ttotal_regions\\n\"; } "
            f"else {{ print \"$_\\t{target_total}\\n\"; }}'"
        )
    else:
        usage()

    subprocess.run(giggle_cmd, shell=True, executable='/bin/bash')

    # Remove temporary files
    os.remove(tmp_bed_name)
    os.remove(tmp_bed_gz)

if __name__ == "__main__":
    main()