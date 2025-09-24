#!/usr/bin/env python3

import sys
import csv
import gzip
import os

def parse_info(info_field):
    """Parse INFO field into a dictionary"""
    info_dict = {}
    for entry in info_field.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[entry] = True
    return info_dict

def parse_vcf(vcf_path, tsv_path, bed_path=None):
    info_keys = set()
    records = []
    bed_records = []

    with gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            info_data = parse_info(fields[7])
            info_keys.update(info_data.keys())

            record = {
                'CHROM': fields[0],
                'POS': fields[1],
                'ID': fields[2],
                'REF': fields[3],
                'ALT': fields[4],
                'QUAL': fields[5],
                'FILTER': fields[6],
                **info_data
            }
            records.append(record)

            # Prepare BED fields (optional)
            if bed_path:
                chrom = fields[0]
                start = int(fields[1]) #- 1  # BED is 0-based
                end = int(info_data.get('END', fields[1]))  # fallback to POS if no END
                svtype = info_data.get('SVTYPE', 'NA')
                bed_records.append((chrom, start, end, svtype))

    # Write TSV
    all_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'] + sorted(info_keys)
    with open(tsv_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=all_columns, delimiter='\t')
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    # Write BED (optional)
    if bed_path:
        with open(bed_path, 'w') as bed:
            for chrom, start, end, svtype in bed_records:
                bed.write(f"{chrom}\t{start}\t{end}\t{svtype}\n")

# if __name__ == "__main__":
#     if 'get_ipython' in globals():
#         # Debug/test in Spyder
#         input_vcf = "/media/rj/SSD_500GB/CNV_OVARIA/PLAZMY/variant_calls/3418-23_plazma/cnvkit/CNV_calls.vcf"
#         output_tsv = "/media/rj/SSD_500GB/CNV_OVARIA/PLAZMY/variant_calls/3418-23_plazma/cnvkit/CNV_calls.test.tsv"
#         output_bed = "/media/rj/SSD_500GB/CNV_OVARIA/PLAZMY/variant_calls/3418-23_plazma/cnvkit/CNV_calls.test.bed"
#         parse_vcf(input_vcf, output_tsv, output_bed)
#     else:
#         if len(sys.argv) < 3:
#             print(f"Usage: {sys.argv[0]} input.vcf[.gz] output.tsv [output.bed]", file=sys.stderr)
#             sys.exit(1)

#         input_vcf = sys.argv[1]
#         output_tsv = sys.argv[2]
#         output_bed = sys.argv[3] if len(sys.argv) > 3 else None

#         parse_vcf(input_vcf, output_tsv, output_bed)

        
if __name__ == "__main__":
        if len(sys.argv) < 3:
            print(f"Usage: {sys.argv[0]} input.vcf[.gz] output.tsv [output.bed]", file=sys.stderr)
            sys.exit(1)

        input_vcf = sys.argv[1]
        output_tsv = sys.argv[2]
        output_bed = sys.argv[3] if len(sys.argv) > 3 else None

        parse_vcf(input_vcf, output_tsv, output_bed)
