#!/usr/bin/env python3
import sys, csv, gzip, os, collections

def parse_info(info_field):
    info_dict = {}
    if not info_field:
        return info_dict
    for entry in info_field.split(';'):
        entry = entry.strip()
        if not entry:
            continue
        if '=' in entry:
            key, value = entry.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[entry] = True
    return info_dict

def parse_vcf(vcf_path, tsv_path, bed_path=None, debug=True):
    info_keys = set()
    records = []
    bed_records = []
    chrom_counter = collections.Counter()

    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt', newline='') as vcf:
        for raw in vcf:
            line = raw.strip()
            if not line:
                continue
            if line.lstrip().startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 8:
                fields = line.split()

            if len(fields) < 8:
                if debug:
                    sys.stderr.write(f"[WARN] Line with < 8 fields: {line[:80]}...\n")
                continue

            chrom = fields[0].strip()
            pos   = fields[1].strip()
            info  = fields[7].strip()

            chrom_counter[chrom] += 1

            info_data = parse_info(info)
            info_keys.update(info_data.keys())

            record = {
                'CHROM': chrom,
                'POS': pos,
                'ID': fields[2],
                'REF': fields[3],
                'ALT': fields[4],
                'QUAL': fields[5],
                'FILTER': fields[6],
                **info_data
            }
            records.append(record)

            if bed_path:
                try:
                    start0 = max(0, int(pos) - 1)
                except ValueError:
                    start0 = 0
                try:
                    end = int(info_data.get('END', pos))
                except ValueError:
                    end = start0 + 1
                svtype = info_data.get('SVTYPE', 'NA')
                bed_records.append((chrom, start0, end, svtype))

    # escribir TSV
    all_columns = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER'] + sorted(info_keys)
    with open(tsv_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=all_columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for rec in records:
            for k in all_columns:
                rec.setdefault(k, '')
            writer.writerow(rec)

    if bed_path:
        with open(bed_path, 'w') as bed:
            for chrom, start0, end, svtype in bed_records:
                bed.write(f"{chrom}\t{start0}\t{end}\t{svtype}\n")

    if debug:
        sys.stderr.write("[INFO] Chromosomes detected and counted:\n")
        for chrom, cnt in chrom_counter.most_common():
            sys.stderr.write(f"  {chrom}: {cnt}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} input.vcf[.gz] output.tsv [output.bed]", file=sys.stderr)
        sys.exit(1)
    parse_vcf(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv)>3 else None)