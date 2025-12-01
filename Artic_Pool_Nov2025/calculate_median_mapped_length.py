#!/usr/bin/env python3
import sys
import re
import statistics

def parse_cigar_length(cigar):
    """Calculate mapped length from CIGAR (M, =, X operations only)"""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    length = sum(int(num) for num, op in ops if op in ['M', '=', 'X'])
    return length

lengths = []

for line in sys.stdin:
    if line.startswith('@'):
        continue
    fields = line.strip().split('\t')
    if len(fields) < 6:
        continue

    flag = int(fields[1])
    cigar = fields[5]

    # Skip unmapped reads
    if flag & 0x4:
        continue

    length = parse_cigar_length(cigar)
    if length > 0:
        lengths.append(length)

if lengths:
    median = statistics.median(lengths)
    mean = statistics.mean(lengths)
    print(f"Reads analyzed: {len(lengths):,}")
    print(f"Median mapped length: {median:.1f} bp")
    print(f"Mean mapped length: {mean:.1f} bp")
    print(f"Min: {min(lengths)} bp")
    print(f"Max: {max(lengths)} bp")
else:
    print("No mapped reads found")
