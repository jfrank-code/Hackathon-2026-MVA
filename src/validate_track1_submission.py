#!/usr/bin/env python3
"""Structural validator for the MVA Track 1 CSV. Does not validate biology or score."""
import argparse, csv, math, sys
from pathlib import Path
REQUIRED = ['proband_id','chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2','epcr','finding_type','notes']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('csv_path')
    args=ap.parse_args()
    p=Path(args.csv_path)
    if not p.exists(): raise SystemExit(f'ERROR: file not found: {p}')
    with p.open(encoding='utf-8-sig', newline='') as f:
        reader=csv.DictReader(f); rows=list(reader); cols=reader.fieldnames or []
    errors=[]; warnings=[]
    missing=[c for c in REQUIRED if c not in cols]
    if missing: errors.append('missing required columns: '+', '.join(missing))
    if not 1 <= len(rows) <= 10: errors.append(f'rows must be 1..10, observed {len(rows)}')
    seen=set(); epcrs=[]
    for i,r in enumerate(rows,1):
        if r.get('proband_id','').strip()!='PROBAND01': errors.append(f'row {i}: proband_id must be PROBAND01')
        try:
            epcr=float(r.get('epcr',''))
            if not math.isfinite(epcr) or not (0 < epcr <= 1): errors.append(f'row {i}: epcr not in (0,1]')
            epcrs.append(epcr)
        except ValueError: errors.append(f'row {i}: EPCR is not numeric')
        if r.get('finding_type','').strip().lower() not in {'primary','secondary'}: errors.append(f'row {i}: finding_type must be primary/secondary')
        chrom1=r.get('chrom_1','').strip(); chrom2=r.get('chrom_2','').strip()
        if not chrom1.startswith('chr'): errors.append(f'row {i}: chrom_1 must have chr prefix')
        for idx in (1,2):
            chrom=r.get(f'chrom_{idx}','').strip(); pos=r.get(f'pos_{idx}','').strip(); ref=r.get(f'ref_{idx}','').strip(); alt=r.get(f'alt_{idx}','').strip()
            if idx==2 and not any((chrom,pos,ref,alt)): continue
            if not all((chrom,pos,ref,alt)): errors.append(f'row {i}: variant {idx} is partially populated')
            if chrom and not chrom.startswith('chr'): errors.append(f'row {i}: chrom_{idx} must have chr prefix')
            if pos:
                try:
                    if int(pos)<=0: errors.append(f'row {i}: pos_{idx} must be positive')
                except ValueError: errors.append(f'row {i}: pos_{idx} must be integer')
            if ref and alt and ref.upper()==alt.upper(): errors.append(f'row {i}: REF equals ALT for variant {idx}')
        key=tuple(r.get(k,'').strip().upper() for k in ['chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2'])
        if key in seen: errors.append(f'row {i}: duplicate finding')
        seen.add(key)
    if len(epcrs)>1 and any(epcrs[i] < epcrs[i+1] for i in range(len(epcrs)-1)):
        warnings.append('EPCR not in descending order; official evaluator may sort, but a ranked file is clearer')
    if any('not calibrated' in r.get('notes','').lower() for r in rows): warnings.append('notes disclose EPCR limitation; verify notes allowed by current form')
    print('File:',p); print('Rows:',len(rows)); print('Columns:',','.join(cols))
    if errors:
        print('\nFAIL')
        for x in errors: print('ERROR:',x)
        return 1
    print('\nPASS structural checks only — biology, exact evaluator score, GRCh38 and EPCR calibration are NOT validated here.')
    for x in warnings: print('WARNING:',x)
    return 0
if __name__=='__main__': sys.exit(main())
