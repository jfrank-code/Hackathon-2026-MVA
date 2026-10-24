#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VariantForge v0.5.4 — association-aware re-ranking of PAIR-TRACE output.

This is a post-processing layer for the existing v0.5.3 pair table. It does not
re-open/re-process the VCF. It preserves source columns, annotates the disease-
association class, re-ranks candidate pairs, flags nearby indels for review, and
writes a capped Track 1 draft whose EPCR values are explicitly placeholders.

It does NOT classify variants, prove phase, infer diagnosis, or produce calibrated
causal probabilities. Review all outputs before use.
"""
from __future__ import annotations
import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Any

DEFAULT_BASE = Path('/mnt/c/MVA_Hackathon/processed')
DEFAULT_INPUT = DEFAULT_BASE / 'variantforge_v0_5_3' / 'variantforge_v053_pair_candidates.tsv'
DEFAULT_OUT = DEFAULT_BASE / 'variantforge_v0_5_4'
MISSING = {'', 'NA', 'N/A', '.', 'NONE', 'NULL'}


def clean(x: Any) -> str:
    if x is None:
        return ''
    s = str(x).strip()
    return '' if s.upper() in MISSING else s


def number(x: Any, default: float = 0.0) -> float:
    try:
        z = float(clean(x))
        return z if math.isfinite(z) else default
    except (TypeError, ValueError):
        return default


def integer(x: Any, default: int = 0) -> int:
    return int(number(x, default))


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding='utf-8-sig', newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        rows = list(r)
        return rows, list(r.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, '') if row.get(c) is not None else 'NA' for c in columns})


def norm_disease_id(value: str) -> str:
    # Accept disease IDs embedded among semicolon/pipe-separated identifiers.
    s = clean(value).upper().replace(' ', '')
    return s


def create_default_override(path: Path) -> None:
    """Write transparent, source-backed classification for OMIM:212750.

    OMIM:212750 is CELIAC DISEASE, SUSCEPTIBILITY TO, 1. It is a complex
    susceptibility association, not evidence that two variants form a Mendelian
    causal pair. This only downweights the phenotype association for pair ranking;
    it does not remove HLA variants or declare them benign.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['DISEASE_ID', 'ASSOCIATION_CLASS', 'PENALTY', 'SOURCE', 'REVIEWED'])
        w.writerow([
            'OMIM:212750', 'COMPLEX_SUSCEPTIBILITY', '1',
            'GeneReviews Celiac Disease, Table B and Molecular Pathogenesis: https://www.ncbi.nlm.nih.gov/books/NBK1727/',
            'YES'
        ])


def load_overrides(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        create_default_override(path)
    result: dict[str, dict[str, str]] = {}
    with path.open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            disease = norm_disease_id(row.get('DISEASE_ID', ''))
            if not disease:
                continue
            reviewed = clean(row.get('REVIEWED', '')).upper() in {'YES', 'Y', 'TRUE', '1'}
            if not reviewed:
                continue
            result[disease] = {
                'class': clean(row.get('ASSOCIATION_CLASS', 'UNCLASSIFIED')) or 'UNCLASSIFIED',
                'penalty': str(max(0, min(5, integer(row.get('PENALTY'), 0)))),
                'source': clean(row.get('SOURCE', '')),
            }
    return result


def classify_association(disease_id: str, overrides: dict[str, dict[str, str]]) -> dict[str, Any]:
    key = norm_disease_id(disease_id)
    # Use exact id matching, never substring matching.
    if key in overrides:
        ov = overrides[key]
        return {'class': ov['class'], 'penalty': integer(ov['penalty']), 'source': ov['source'], 'reviewed': 'YES'}
    if not key:
        return {'class': 'NO_DISEASE_ASSOCIATION_ID', 'penalty': 0, 'source': '', 'reviewed': 'NO'}
    return {'class': 'UNCLASSIFIED_REQUIRES_REVIEW', 'penalty': 0, 'source': '', 'reviewed': 'NO'}


def variant_region(chrom: str, pos: str) -> str:
    """Broad chr6 MHC cluster for submission diversification only (GRCh38)."""
    c = clean(chrom).lower()
    p = integer(pos, -1)
    if c in {'chr6', '6'} and 25_000_000 <= p <= 34_000_000:
        return 'MHC_REGION_CHR6_25_34MB'
    return f'{c}:{p // 5_000_000}MB_BIN' if p >= 0 and c else 'UNKNOWN_REGION'


def nearby_indel_review(row: dict[str, Any]) -> tuple[str, str]:
    c1, c2 = clean(row.get('CHROM_1')), clean(row.get('CHROM_2'))
    p1, p2 = integer(row.get('POS_1'), -1), integer(row.get('POS_2'), -1)
    r1, a1 = clean(row.get('REF_1')), clean(row.get('ALT_1'))
    r2, a2 = clean(row.get('REF_2')), clean(row.get('ALT_2'))
    if p1 < 0 or p2 < 0 or not c1 or c1 != c2:
        return 'NO_SAME_CHROMOSOME_DISTANCE', 'NA'
    e1 = p1 + max(1, len(r1)) - 1
    e2 = p2 + max(1, len(r2)) - 1
    gap = max(0, max(p1, p2) - min(e1, e2) - 1)
    is_indel = len(r1) != len(a1) or len(r2) != len(a2)
    if is_indel and gap <= 1000:
        return 'REVIEW_NEARBY_INDEL_NOT_PROVEN_COMPLEX_EVENT', str(gap)
    return 'NO_NEARBY_INDEL_FLAG', str(gap)


def rank_tuple(row: dict[str, Any]) -> tuple[Any, ...]:
    # Evidence-first ordinal ordering, not a calibrated causal score.
    # For a reviewed complex-susceptibility association, HPO counts from that
    # association are retained raw but are not used as Mendelian pair evidence.
    assoc_class = row.get('HPO_ASSOCIATION_CLASS', '')
    hpo_suppressed = assoc_class == 'COMPLEX_SUSCEPTIBILITY'
    hpo_subtype = 0 if hpo_suppressed else integer(row.get('HPO_SPECIFIC_SUBTYPE_MATCHES'))
    hpo_exact = 0 if hpo_suppressed else integer(row.get('HPO_EXACT_MATCHES'))
    hpo_relation = 0 if hpo_suppressed else integer(row.get('HPO_ANY_RELATION_COUNT'))
    conflicted = sum(clean(row.get(c)).upper() == 'CONFLICTING' for c in ('CLINVAR_STATE_1', 'CLINVAR_STATE_2'))
    nearby_flag = row.get('NEARBY_INDEL_REVIEW_FLAG') == 'REVIEW_NEARBY_INDEL_NOT_PROVEN_COMPLEX_EVENT'
    return (
        -integer(row.get('MODEL_COMPATIBILITY_AXIS')),
        integer(row.get('CLINVAR_BENIGN_COUNT')),
        conflicted,
        -integer(row.get('CLINVAR_PLP_COUNT')),
        integer(row.get('ASSOCIATION_CLASS_PENALTY')),
        -hpo_subtype,
        -hpo_exact,
        -hpo_relation,
        -integer(row.get('LOF_VARIANT_COUNT')),
        -integer(row.get('RARE_AF_COUNT')),
        -integer(row.get('TECHNICAL_PASS_COUNT')),
        -integer(row.get('GNOMAD_AF_KNOWN_COUNT')),
        1 if nearby_flag else 0,
        clean(row.get('GENE')),
        clean(row.get('VARIANT_1_KEY')),
        clean(row.get('VARIANT_2_KEY')),
    )


def submission_region(row: dict[str, Any]) -> str:
    # If either variant is in the chr6 MHC span, group the pair as one region.
    r1 = variant_region(clean(row.get('CHROM_1')), clean(row.get('POS_1')))
    r2 = variant_region(clean(row.get('CHROM_2')), clean(row.get('POS_2')))
    if 'MHC_REGION_CHR6_25_34MB' in {r1, r2}:
        return 'MHC_REGION_CHR6_25_34MB'
    # The pair is intragenic in this table; use the first locus bin as a conservative cluster.
    return r1


def main() -> int:
    ap = argparse.ArgumentParser(description='VariantForge v0.5.4 association-aware re-ranking; no VCF reprocessing.')
    ap.add_argument('--input', default=str(DEFAULT_INPUT), help='v0.5.3 pair_candidates TSV')
    ap.add_argument('--out', default=str(DEFAULT_OUT), help='Output folder')
    ap.add_argument('--overrides', default=None, help='Optional reviewed disease-association override TSV')
    ap.add_argument('--submission-candidates', type=int, default=10, help='Rows in review draft (1-10)')
    ap.add_argument('--max-pairs-per-gene', type=int, default=1, help='Max pairs per gene in draft')
    ap.add_argument('--draft-epcr', type=float, default=0.95, help='ORDER PLACEHOLDER only, not a probability estimate')
    args = ap.parse_args()
    if not 1 <= args.submission_candidates <= 10:
        ap.error('--submission-candidates must be between 1 and 10')
    if args.max_pairs_per_gene < 1:
        ap.error('--max-pairs-per-gene must be >= 1')
    if not 0 < args.draft_epcr <= 1 or args.draft_epcr - 0.05*(args.submission_candidates-1) <= 0:
        ap.error('--draft-epcr must leave all descending placeholders in (0,1]')

    in_path = Path(args.input).expanduser()
    out = Path(args.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    if not in_path.exists():
        raise FileNotFoundError(f'Input pair table not found: {in_path}')
    override_path = Path(args.overrides).expanduser() if args.overrides else out / 'disease_association_overrides.tsv'
    overrides = load_overrides(override_path)
    rows, columns = read_tsv(in_path)
    if not rows:
        raise RuntimeError(f'Input pair table has no rows: {in_path}')

    added_columns = [
        'HPO_ASSOCIATION_CLASS', 'ASSOCIATION_CLASS_PENALTY', 'HPO_EFFECTIVE_EXACT_FOR_PAIR_RANK',
        'HPO_EFFECTIVE_SUBTYPE_FOR_PAIR_RANK', 'HPO_ASSOCIATION_OVERRIDE_REVIEWED',
        'HPO_ASSOCIATION_OVERRIDE_SOURCE', 'NEARBY_INDEL_REVIEW_FLAG', 'PAIR_GAP_BP',
        'MHC_REGION_CLUSTER',
    ]
    work: list[dict[str, Any]] = []
    for original in rows:
        r: dict[str, Any] = dict(original)
        assoc = classify_association(clean(r.get('HPO_BEST_DISEASE_ID')), overrides)
        r['HPO_ASSOCIATION_CLASS'] = assoc['class']
        r['ASSOCIATION_CLASS_PENALTY'] = assoc['penalty']
        r['HPO_ASSOCIATION_OVERRIDE_REVIEWED'] = assoc['reviewed']
        r['HPO_ASSOCIATION_OVERRIDE_SOURCE'] = assoc['source']
        if assoc['class'] == 'COMPLEX_SUSCEPTIBILITY':
            r['HPO_EFFECTIVE_EXACT_FOR_PAIR_RANK'] = 0
            r['HPO_EFFECTIVE_SUBTYPE_FOR_PAIR_RANK'] = 0
            r['AUDIT_FLAGS'] = ';'.join(x for x in [clean(r.get('AUDIT_FLAGS')), 'HPO_ASSOCIATION_COMPLEX_SUSCEPTIBILITY_NOT_MENDELIAN_PAIR_EVIDENCE'] if x)
        else:
            r['HPO_EFFECTIVE_EXACT_FOR_PAIR_RANK'] = integer(r.get('HPO_EXACT_MATCHES'))
            r['HPO_EFFECTIVE_SUBTYPE_FOR_PAIR_RANK'] = integer(r.get('HPO_SPECIFIC_SUBTYPE_MATCHES'))
        flag, gap = nearby_indel_review(r)
        r['NEARBY_INDEL_REVIEW_FLAG'] = flag
        r['PAIR_GAP_BP'] = gap
        r['MHC_REGION_CLUSTER'] = 'YES' if submission_region(r) == 'MHC_REGION_CHR6_25_34MB' else 'NO'
        work.append(r)

    work.sort(key=rank_tuple)
    for rank, r in enumerate(work, 1):
        r['FINAL_PRIORITY_RANK_V054'] = rank
        r['RANKING_POLICY_V054'] = 'ASSOCIATION_AWARE_EVIDENCE_FIRST_ORDINAL_NOT_CALIBRATED'
    output_columns = columns + [c for c in added_columns if c not in columns] + [
        c for c in ['FINAL_PRIORITY_RANK_V054', 'RANKING_POLICY_V054'] if c not in columns
    ]
    pair_out = out / 'variantforge_v054_pair_candidates.tsv'
    write_tsv(pair_out, work, output_columns)

    # Diverse review draft: at most one pair per gene and one per broad MHC region.
    chosen: list[dict[str, Any]] = []
    gene_counts: Counter[str] = Counter()
    used_regions: set[str] = set()
    rejected_mhc_diversity = 0
    for r in work:
        gene = clean(r.get('GENE'))
        region = submission_region(r)
        is_mhc = region == 'MHC_REGION_CHR6_25_34MB'
        if gene_counts[gene] >= args.max_pairs_per_gene:
            continue
        if is_mhc and region in used_regions:
            rejected_mhc_diversity += 1
            continue
        chosen.append(r)
        gene_counts[gene] += 1
        if is_mhc:
            used_regions.add(region)
        if len(chosen) >= args.submission_candidates:
            break

    sub_path = out / 'variantforge_v054_track1_submission_draft.csv'
    with sub_path.open('w', encoding='utf-8', newline='') as f:
        cols = ['proband_id','chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2','epcr','finding_type','notes']
        w = csv.DictWriter(f, fieldnames=cols, lineterminator='\n')
        w.writeheader()
        for i, r in enumerate(chosen, 1):
            w.writerow({
                'proband_id': 'PROBAND01',
                'chrom_1': r.get('CHROM_1',''), 'pos_1': r.get('POS_1',''), 'ref_1': r.get('REF_1',''), 'alt_1': r.get('ALT_1',''),
                'chrom_2': r.get('CHROM_2',''), 'pos_2': r.get('POS_2',''), 'ref_2': r.get('REF_2',''), 'alt_2': r.get('ALT_2',''),
                'epcr': f'{args.draft_epcr - 0.05*(i-1):.2f}', 'finding_type': 'primary',
                'notes': f'v0.5.4 review-order placeholder #{i}; association-aware ranking; EPCR NOT CALIBRATED; validate before submission'
            })

    report_path = out / 'RUN_REPORT.txt'
    with report_path.open('w', encoding='utf-8') as f:
        f.write('VariantForge v0.5.4 — association-aware re-ranking\n')
        f.write('='*68 + '\n')
        f.write(f'Input pairs: {in_path}\nOutput: {out}\nReviewed overrides: {override_path}\nRows: {len(work)}\n')
        f.write('This layer reads the v0.5.3 pair TSV only; it does not reprocess VCF data.\n')
        f.write('OMIM:212750 is annotated as COMPLEX_SUSCEPTIBILITY based on GeneReviews Celiac Disease (Table B and Molecular Pathogenesis): https://www.ncbi.nlm.nih.gov/books/NBK1727/\n')
        f.write('This annotation suppresses that association’s HPO counts as support for a Mendelian pair; it does NOT delete HLA variants or classify them as benign.\n')
        f.write('Nearby indel flag is an audit flag only; proximity does not prove a complex variant event.\n')
        f.write('MHC-region limit in the draft is a diversification heuristic only; candidates remain in the full ranking.\n')
        f.write('The EPCR draft values are placeholders for order and are NOT calibrated causal probabilities. Do not submit them as probability estimates.\n')
        f.write(f'Draft rows: {len(chosen)}; one pair per gene and at most one pair from chr6:25-34 Mb GRCh38; candidates skipped by MHC diversity: {rejected_mhc_diversity}\n')
        f.write('\nTOP 30 PAIRS\n')
        for r in work[:30]:
            f.write(
                f"{r['FINAL_PRIORITY_RANK_V054']} | {r.get('GENE')} | disease={r.get('HPO_BEST_DISEASE_ID')} | class={r.get('HPO_ASSOCIATION_CLASS')} | model={r.get('INHERITANCE_MODEL_HINT')} | reviewed={r.get('INHERITANCE_HINT_REVIEWED')} | PLP={r.get('CLINVAR_PLP_COUNT')} | benign={r.get('CLINVAR_BENIGN_COUNT')} | LoF={r.get('LOF_VARIANT_COUNT')} | HPO raw={r.get('HPO_EXACT_MATCHES')} effective={r.get('HPO_EFFECTIVE_EXACT_FOR_PAIR_RANK')} | nearby indel={r.get('NEARBY_INDEL_REVIEW_FLAG')} | {r.get('VARIANT_1_KEY')} + {r.get('VARIANT_2_KEY')}\n"
            )
        f.write('\nDRAFT SUBMISSION ROWS\n')
        for i, r in enumerate(chosen, 1):
            f.write(f"{i}. {r.get('GENE')} | {r.get('VARIANT_1_KEY')} + {r.get('VARIANT_2_KEY')}\n")

    print('VariantForge v0.5.4 — asociación-aware re-ranking completado')
    print(f'Pares leídos y reclasificados: {len(work)}')
    print(f'Salida: {out}')
    print('Archivos: variantforge_v054_pair_candidates.tsv, variantforge_v054_track1_submission_draft.csv, RUN_REPORT.txt, disease_association_overrides.tsv')
    print('\nTOP 15')
    for r in work[:15]:
        print(
            f"{r['FINAL_PRIORITY_RANK_V054']:>4} | {clean(r.get('GENE')):<12} | "
            f"Assoc={clean(r.get('HPO_ASSOCIATION_CLASS')):<24} | "
            f"Modelo={clean(r.get('INHERITANCE_MODEL_HINT')):<18} | "
            f"P/LP={integer(r.get('CLINVAR_PLP_COUNT'))} | LoF={integer(r.get('LOF_VARIANT_COUNT'))} | "
            f"HPO raw/effective={integer(r.get('HPO_EXACT_MATCHES'))}/{integer(r.get('HPO_EFFECTIVE_EXACT_FOR_PAIR_RANK'))} | "
            f"{clean(r.get('VARIANT_1_KEY'))} + {clean(r.get('VARIANT_2_KEY'))}"
        )
    print('\nBORRADOR CSV — revisar; EPCR son placeholders, no probabilidades calibradas')
    for i, r in enumerate(chosen, 1):
        print(f"{i}. {r.get('GENE')} | {r.get('VARIANT_1_KEY')} + {r.get('VARIANT_2_KEY')}")
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, csv.Error, ValueError, RuntimeError) as e:
        raise SystemExit(f'ERROR: {e}')
