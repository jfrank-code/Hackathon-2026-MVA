#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PIVOT-MVA v1.0 — Pairwise Inheritance + Phenotype-domain + Traceable evidence.

Lightweight post-processor for existing VariantForge TSVs. It does not read or
reprocess the raw VCF. Python standard library only.

Important: this is an ordinal evidence-ranking model, not a diagnostic system
and not an empirically calibrated probability model. EPCR fields in the output
are explicit submission estimates/placeholders and must be reviewed.
"""
from __future__ import annotations
import argparse, csv, math, re, sys
from collections import defaultdict, Counter
from pathlib import Path
from typing import Any

MISSING = {'', 'NA', 'N/A', '.', 'NONE', 'NULL', 'NAN'}
DEFAULT_ROOT = Path('/mnt/c/MVA_Hackathon/processed')
DEFAULT_PAIRS = DEFAULT_ROOT / 'variantforge_v0_5_3' / 'variantforge_v053_pair_candidates.tsv'
DEFAULT_HPO = DEFAULT_ROOT / 'variantforge_v0_4_1' / 'gene_disease_summary.tsv'
DEFAULT_INHERITANCE = DEFAULT_ROOT / 'variantforge_v0_5' / 'inheritance_overrides_curated.tsv'
DEFAULT_OUT = DEFAULT_ROOT / 'pivot_mva_v1_0'

# IC values taken from the user's already-generated HPO report; scores are
# relative phenotype-support descriptors, never probabilities.
PHENOTYPES = {
    'rhabdomyosarcoma': ('oncologic', 6.2328),
    'nephrocalcinosis': ('renal', 4.5080),
    'short stature': ('growth', 1.6942),
    'failure to thrive': ('growth', 2.1723),
    'skeletal muscle atrophy': ('neuromuscular', 2.5125),
    'premature birth': ('perinatal', 4.1259),
    'small for gestational age': ('perinatal', 3.5359),
}
DOMAIN_ORDER = ('oncologic', 'renal', 'growth', 'neuromuscular', 'perinatal')

# Disease-level classifications only. The algorithm has no hard-coded variant
# coordinates. Unknown associations remain unknown; they are not assumed to be
# monogenic. BUB1B/MVA1 is a curated Mendelian relationship; celiac susceptibility
# is explicitly complex. Sources are stored with the configuration and report.
DEFAULT_ASSOCIATIONS = [
    {
        'DISEASE_ID': 'OMIM:257300', 'ASSOCIATION_CLASS': 'MENDELIAN_DISEASE',
        'INHERITANCE': 'AR', 'REVIEWED': 'YES',
        'SOURCE': 'ClinGen gene-disease validity for BUB1B/MVA1: https://search.clinicalgenome.org/kb/gene-validity/CGGV%3Aassertion_59147f27-d5a3-4760-ba8d-0429bae3c906-2019-11-22T14%3A53%3A26.352Z',
        'NOTE': 'Disease-level rule; pair phase remains unknown.'
    },
    {
        'DISEASE_ID': 'OMIM:212750', 'ASSOCIATION_CLASS': 'COMPLEX_SUSCEPTIBILITY',
        'INHERITANCE': 'NOT_APPLICABLE', 'REVIEWED': 'YES',
        'SOURCE': 'GeneReviews Celiac Disease: https://www.ncbi.nlm.nih.gov/books/NBK1727/',
        'NOTE': 'Do not use this association as Mendelian-pair support.'
    },
]


def clean(v: Any) -> str:
    if v is None: return ''
    s = str(v).strip()
    return '' if s.upper() in MISSING else s


def num(v: Any, default: float | None = None) -> float | None:
    s = clean(v)
    if not s: return default
    try:
        x = float(s)
        return x if math.isfinite(x) else default
    except (ValueError, TypeError): return default


def integer(v: Any, default: int = 0) -> int:
    x = num(v, None)
    return int(x) if x is not None else default


def truth(v: Any) -> bool:
    return clean(v).upper() in {'YES', 'Y', 'TRUE', '1', 'PASS'}


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding='utf-8-sig', newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        rows = list(r)
        return rows, list(r.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    cols: list[str] = []
    for r in rows:
        for c in r:
            if c not in cols: cols.append(c)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({c: ('NA' if r.get(c) is None else r.get(c, '')) for c in cols})


def norm_id(s: str) -> str:
    return re.sub(r'\s+', '', clean(s)).upper()


def norm_gene(s: str) -> str:
    return clean(s).upper()


def make_defaults(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        cols = ['DISEASE_ID','ASSOCIATION_CLASS','INHERITANCE','REVIEWED','SOURCE','NOTE']
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(DEFAULT_ASSOCIATIONS)


def load_associations(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists(): make_defaults(path)
    out = {}
    with path.open(encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            did = norm_id(r.get('DISEASE_ID', ''))
            if did and truth(r.get('REVIEWED', '')):
                out[did] = {k: clean(v) for k, v in r.items()}
    return out


def load_inheritance(path: Path | None) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not path or not path.exists(): return out
    with path.open(encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            gene = norm_gene(r.get('GENE', ''))
            if gene and truth(r.get('REVIEWED', '')):
                out[gene] = {k: clean(v) for k, v in r.items()}
    return out


def hpo_index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    out = {}
    for r in rows:
        gene = norm_gene(r.get('GENE', ''))
        did = norm_id(r.get('DISEASE_ID', ''))
        if gene and did: out[(gene, did)] = r
    return out


def phenotype_parts(raw: str) -> list[tuple[str, str]]:
    """Return (term, relation) from 'term=EXACT; term2=MORE_SPECIFIC_SUBTYPE'."""
    out = []
    for part in clean(raw).split(';'):
        if '=' not in part: continue
        term, rel = part.rsplit('=', 1)
        term = re.sub(r'\s+', ' ', term.strip()).lower()
        rel = clean(rel).upper()
        if term: out.append((term, rel))
    return out


def parse_hpo_profile(row: dict[str, str] | None, assoc_class: str) -> tuple[dict[str, float], dict[str, str], int, int, str]:
    components = {d: 0.0 for d in DOMAIN_ORDER}
    matched = {}
    if not row or assoc_class == 'COMPLEX_SUSCEPTIBILITY':
        return components, matched, 0, 0, ''
    for term, relation in phenotype_parts(row.get('MATCHED_PHENOTYPES', '')):
        if term not in PHENOTYPES: continue
        domain, ic = PHENOTYPES[term]
        if relation == 'EXACT': factor = 1.0
        elif relation in {'MORE_SPECIFIC_SUBTYPE', 'SPECIFIC_SUBTYPE', 'SUBTYPE'}: factor = 0.65
        elif relation in {'BROAD_ANCESTOR_ONLY', 'ANCESTOR', 'BROAD_ANCESTOR'}: factor = 0.25
        else: factor = 0.4
        value = ic * factor
        # Collapse correlated terms: within each phenotype domain, only the
        # strongest informative match contributes to the aggregate.
        if value > components[domain]:
            components[domain] = value
        matched[term] = relation
    n_exact = sum(1 for _, rel in phenotype_parts(row.get('MATCHED_PHENOTYPES', '')) if rel == 'EXACT')
    n_rel = sum(1 for _, rel in phenotype_parts(row.get('MATCHED_PHENOTYPES', '')) if rel)
    return components, matched, n_exact, n_rel, clean(row.get('MATCHED_PHENOTYPES', ''))


def genotype_class(gt: str) -> str:
    g = clean(gt).replace('|', '/').upper()
    if g in {'0/1', '1/0'}: return 'HET'
    if g in {'1/1'}: return 'HOM_ALT'
    if g in {'0/0'}: return 'REF'
    return 'UNKNOWN'


def variant_af(row: dict[str, str], idx: int) -> float | None:
    # No exact AF remains missing, not zero.
    for c in (f'GNOMAD_AF_{idx}', f'GNOMAD_AF_EXACT_{idx}', f'AF_{idx}'):
        if c in row:
            x = num(row.get(c), None)
            if x is not None and 0 <= x <= 1: return x
    return None


def exact_clinvar_state(row: dict[str, str], idx: int) -> str:
    state = clean(row.get(f'CLINVAR_STATE_{idx}', '')).upper()
    exact = clean(row.get(f'CLINVAR_EXACT_{idx}', '')).upper()
    # If explicit exact-match is available, require it. Older output may encode
    # exact P/LP count aggregate; fall back to the aggregate only if no per-allele state exists.
    if exact and exact not in {'YES', 'Y', 'TRUE', '1'}:
        return 'NO_EXACT_ALLELE_MATCH_OR_UNCONFIRMED'
    return state


def is_plp(state: str) -> bool:
    return state in {'PATHOGENIC', 'LIKELY_PATHOGENIC', 'PATHOGENIC_OR_LIKELY_PATHOGENIC', 'PATHOGENIC/LIKELY_PATHOGENIC'}


def build_features(r: dict[str, str], hpomap: dict[tuple[str, str], dict[str, str]], associations: dict[str, dict[str, str]], inheritance: dict[str, dict[str, str]]) -> dict[str, Any]:
    row = dict(r)
    gene = norm_gene(row.get('GENE', ''))
    did = norm_id(row.get('HPO_BEST_DISEASE_ID', ''))
    assoc = associations.get(did, {})
    assoc_class = clean(assoc.get('ASSOCIATION_CLASS', 'UNCLASSIFIED_REQUIRES_REVIEW')) or 'UNCLASSIFIED_REQUIRES_REVIEW'
    assoc_reviewed = truth(assoc.get('REVIEWED', ''))
    source = clean(assoc.get('SOURCE', ''))
    hpo_row = hpomap.get((gene, did))
    if hpo_row is None and not did:
        # Fallback is allowed only when the input has no disease ID. If an ID is
        # present but doesn't join, do not silently borrow a different disease profile.
        same_gene = [v for (g, _), v in hpomap.items() if g == gene]
        if same_gene:
            hpo_row = max(same_gene, key=lambda x: integer(x.get('N_EXACT')) + 0.5*integer(x.get('N_SPECIFIC_SUBTYPE_RELATIONS')))
    domain_scores, matched, parsed_exact, parsed_n, matched_raw = parse_hpo_profile(hpo_row, assoc_class)
    inh_row = inheritance.get(gene, {})
    inh_model = (clean(row.get('INHERITANCE_MODEL_HINT', '')) or clean(inh_row.get('INHERITANCE', '')) or 'UNKNOWN').upper()
    inh_reviewed = truth(row.get('INHERITANCE_HINT_REVIEWED', '')) or truth(inh_row.get('REVIEWED', ''))
    gt1, gt2 = genotype_class(row.get('GT_1', '')), genotype_class(row.get('GT_2', ''))
    if assoc_class == 'COMPLEX_SUSCEPTIBILITY':
        inheritance_tier = 3
        inheritance_status = 'COMPLEX_ASSOCIATION_NOT_MENDELIAN_PAIR_EVIDENCE'
    elif inh_model == 'AR' and inh_reviewed and gt1 == 'HET' and gt2 == 'HET':
        inheritance_tier = 0
        inheritance_status = 'REVIEWED_AR_TWO_HET_PHASE_UNKNOWN'
    elif inh_model == 'AR' and inh_reviewed and 'HOM_ALT' in {gt1, gt2}:
        inheritance_tier = 1
        inheritance_status = 'REVIEWED_AR_GENOTYPES_REQUIRE_REVIEW'
    elif inh_model in {'AD', 'X-LINKED', 'XL', 'X_LINKED'} and inh_reviewed:
        inheritance_tier = 2
        inheritance_status = 'REVIEWED_NON_AR_PAIR_REQUIRES_CONDITION_SPECIFIC_REVIEW'
    else:
        inheritance_tier = 2
        inheritance_status = 'INHERITANCE_UNKNOWN_OR_UNREVIEWED'

    states = [exact_clinvar_state(row, 1), exact_clinvar_state(row, 2)]
    plp_count = sum(is_plp(s) for s in states)
    benign_count = sum(s in {'BENIGN','LIKELY_BENIGN','BENIGN/LIKELY_BENIGN'} for s in states)
    conflicting_count = sum('CONFLICT' in s for s in states)
    vus_count = sum(s in {'VUS','UNCERTAIN_SIGNIFICANCE'} for s in states)
    def norm_condition(value: str) -> str:
        value = clean(value).lower().replace('_', ' ')
        value = re.sub(r'[^a-z0-9 ]+', ' ', value)
        return re.sub(r'\s+', ' ', value).strip()
    condition_candidates = [
        row.get('HPO_BEST_DISEASE_NAME',''), row.get('INHERITANCE_CONDITION_HINT',''),
        inh_row.get('CONDITION',''), assoc.get('CONDITION','')
    ]
    norm_candidates = [norm_condition(x) for x in condition_candidates if clean(x)]
    disease_match_count = 0
    for i, state in enumerate(states, 1):
        if not is_plp(state):
            continue
        disease_labels = [norm_condition(x) for x in re.split(r'[|;]', clean(row.get(f'CLINVAR_DISEASES_{i}', ''))) if clean(x)]
        is_linked = False
        for label in disease_labels:
            for condition in norm_candidates:
                # Require a meaningful label, and match full normalized condition
                # labels rather than gene names or generic terms such as "not provided".
                if len(condition) >= 8 and len(label) >= 8 and (label in condition or condition in label):
                    is_linked = True
        if did == 'OMIM:257300' and any('mosaic variegated aneuploidy' in label for label in disease_labels):
            is_linked = True
        if is_linked:
            disease_match_count += 1

    afs = [variant_af(row, i) for i in (1,2)]
    af_known = sum(a is not None for a in afs)
    rare_count = sum(a is not None and a <= 1e-4 for a in afs)
    nearfix_count = sum(a is not None and a >= 0.95 for a in afs)
    loF_count = integer(row.get('LOF_VARIANT_COUNT'), sum('frameshift' in clean(row.get(f'CONSEQUENCE_{i}','')).lower() or 'stop_gained' in clean(row.get(f'CONSEQUENCE_{i}','')).lower() or 'splice_acceptor' in clean(row.get(f'CONSEQUENCE_{i}','')).lower() or 'splice_donor' in clean(row.get(f'CONSEQUENCE_{i}','')).lower() for i in (1,2)))
    technical_pass = 0
    for i in (1,2):
        dp = num(row.get(f'DP_{i}'), None); gq = num(row.get(f'GQ_{i}'), None)
        afread = num(row.get(f'AF_READS_{i}'), None)
        ok = (dp is None or dp >= 10) and (gq is None or gq >= 20)
        if afread is not None and not 0.20 <= afread <= 0.80: ok = False
        if ok and (dp is not None or gq is not None): technical_pass += 1
    if 'TECHNICAL_PASS_COUNT' in row:
        technical_pass = max(technical_pass, integer(row.get('TECHNICAL_PASS_COUNT')))
    c1, c2 = clean(row.get('CHROM_1','')), clean(row.get('CHROM_2',''))
    p1, p2 = integer(row.get('POS_1'), -1), integer(row.get('POS_2'), -1)
    ref1, alt1 = clean(row.get('REF_1','')), clean(row.get('ALT_1',''))
    ref2, alt2 = clean(row.get('REF_2','')), clean(row.get('ALT_2',''))
    same_chrom = c1.replace('chr','') == c2.replace('chr','') and c1 and c2
    gap = None
    nearby_indel = False
    if same_chrom and p1 >= 0 and p2 >= 0:
        e1 = p1 + max(1,len(ref1))-1; e2 = p2 + max(1,len(ref2))-1
        gap = max(0, max(p1,p2)-min(e1,e2)-1)
        nearby_indel = (len(ref1)!=len(alt1) or len(ref2)!=len(alt2)) and gap <= 1000
    mhc = (c1.lower() in {'chr6','6'} and 25_000_000 <= p1 <= 34_000_000) or (c2.lower() in {'chr6','6'} and 25_000_000 <= p2 <= 34_000_000)
    domains_present = sum(x > 0 for x in domain_scores.values())
    hpo_score = sum(domain_scores.values())
    effective_exact = 0 if assoc_class == 'COMPLEX_SUSCEPTIBILITY' else integer(row.get('HPO_EXACT_MATCHES', parsed_exact))
    effective_subtype = 0 if assoc_class == 'COMPLEX_SUSCEPTIBILITY' else integer(row.get('HPO_SPECIFIC_SUBTYPE_MATCHES'))
    # Clinical linkage of each ClinVar record to target disorder is separate from a global P/LP label.
    plp_unlinked = max(0, plp_count - disease_match_count)
    assoc_tier = 0 if assoc_class == 'MENDELIAN_DISEASE' and assoc_reviewed else (2 if assoc_class == 'COMPLEX_SUSCEPTIBILITY' else 1)
    # `association_class` is disease-level, not a statement that every variant is causal.
    features = {
        **row,
        'PIVOT_ASSOCIATION_CLASS': assoc_class,
        'PIVOT_ASSOCIATION_REVIEWED': 'YES' if assoc_reviewed else 'NO',
        'PIVOT_ASSOCIATION_SOURCE': source,
        'PIVOT_ASSOCIATION_TIER': assoc_tier,
        'PIVOT_INHERITANCE_STATUS': inheritance_status,
        'PIVOT_INHERITANCE_TIER': inheritance_tier,
        'PIVOT_GENOTYPE_1_CLASS': gt1, 'PIVOT_GENOTYPE_2_CLASS': gt2,
        'PIVOT_CLINVAR_PLP_EXACT': plp_count,
        'PIVOT_CLINVAR_DISEASE_MATCH_PLP': disease_match_count,
        'PIVOT_CLINVAR_PLP_UNLINKED': plp_unlinked,
        'PIVOT_CLINVAR_BENIGN_EXACT': benign_count,
        'PIVOT_CLINVAR_CONFLICTING_EXACT': conflicting_count,
        'PIVOT_CLINVAR_VUS_EXACT': vus_count,
        'PIVOT_AF_KNOWN': af_known, 'PIVOT_AF_RARE_LE_1E4': rare_count,
        'PIVOT_AF_NEAR_FIXATION': nearfix_count,
        'PIVOT_LOF_COUNT': loF_count,
        'PIVOT_TECHNICAL_PASS_COUNT': technical_pass,
        'PIVOT_HPO_RAW_EXACT': integer(row.get('HPO_EXACT_MATCHES', parsed_exact)),
        'PIVOT_HPO_EFFECTIVE_EXACT': effective_exact,
        'PIVOT_HPO_EFFECTIVE_SUBTYPE': effective_subtype,
        'PIVOT_HPO_DOMAIN_COUNT': domains_present,
        'PIVOT_HPO_DOMAIN_SCORE': round(hpo_score, 4),
        'PIVOT_HPO_ONCOLOGIC': round(domain_scores['oncologic'],4),
        'PIVOT_HPO_RENAL': round(domain_scores['renal'],4),
        'PIVOT_HPO_GROWTH': round(domain_scores['growth'],4),
        'PIVOT_HPO_NEUROMUSCULAR': round(domain_scores['neuromuscular'],4),
        'PIVOT_HPO_PERINATAL': round(domain_scores['perinatal'],4),
        'PIVOT_HPO_MATCHED_TERMS': matched_raw,
        'PIVOT_HPO_DOMAIN_PROFILE': ';'.join(f'{k}={v:.4f}' for k,v in domain_scores.items() if v > 0) or 'NONE',
        'PIVOT_PHASE_STATUS': clean(row.get('PHASE_STATUS','')) or ('UNKNOWN_NOT_PHASED' if gt1 == 'HET' and gt2 == 'HET' else 'NOT_ASSESSED'),
        'PIVOT_NEARBY_INDEL_REVIEW': 'YES' if nearby_indel else 'NO',
        'PIVOT_PAIR_GAP_BP': gap if gap is not None else 'NA',
        'PIVOT_MHC_REGION_REVIEW': 'YES' if mhc else 'NO',
        'PIVOT_SOURCE_RANK': integer(row.get('FINAL_PRIORITY_RANK_V054'), integer(row.get('FINAL_PRIORITY_RANK'), 999999)),
    }
    return features


def rank_tuple(r: dict[str, Any], omit_domain: str | None = None) -> tuple[Any, ...]:
    # Multi-stage evidence precedence instead of an opaque additive score.
    # Lower tuple is better. Phenotypes are collapsed by domain to limit correlation.
    assoc_tier = integer(r.get('PIVOT_ASSOCIATION_TIER'), 1)
    inh_tier = integer(r.get('PIVOT_INHERITANCE_TIER'), 2)
    domain_scores = {
        'oncologic': num(r.get('PIVOT_HPO_ONCOLOGIC'),0) or 0,
        'renal': num(r.get('PIVOT_HPO_RENAL'),0) or 0,
        'growth': num(r.get('PIVOT_HPO_GROWTH'),0) or 0,
        'neuromuscular': num(r.get('PIVOT_HPO_NEUROMUSCULAR'),0) or 0,
        'perinatal': num(r.get('PIVOT_HPO_PERINATAL'),0) or 0,
    }
    if omit_domain in domain_scores: domain_scores[omit_domain] = 0.0
    hpo_total = sum(domain_scores.values())
    hpo_domain_count = sum(v > 0 for v in domain_scores.values())
    return (
        assoc_tier,
        inh_tier,
        integer(r.get('PIVOT_CLINVAR_BENIGN_EXACT')),
        integer(r.get('PIVOT_CLINVAR_CONFLICTING_EXACT')),
        -integer(r.get('PIVOT_CLINVAR_DISEASE_MATCH_PLP')),
        -integer(r.get('PIVOT_CLINVAR_PLP_EXACT')),
        -round(hpo_total, 4),
        -hpo_domain_count,
        -integer(r.get('PIVOT_HPO_EFFECTIVE_SUBTYPE')),
        -integer(r.get('PIVOT_HPO_EFFECTIVE_EXACT')),
        -integer(r.get('PIVOT_LOF_COUNT')),
        integer(r.get('PIVOT_AF_NEAR_FIXATION')),
        -integer(r.get('PIVOT_AF_RARE_LE_1E4')),
        -integer(r.get('PIVOT_TECHNICAL_PASS_COUNT')),
        1 if r.get('PIVOT_NEARBY_INDEL_REVIEW') == 'YES' else 0,
        1 if r.get('PIVOT_MHC_REGION_REVIEW') == 'YES' else 0,
        clean(r.get('GENE')),
        clean(r.get('VARIANT_1_KEY')),
        clean(r.get('VARIANT_2_KEY')),
    )


def row_pair_key(r: dict[str, Any]) -> tuple[tuple[str,str,str,str], tuple[str,str,str,str]] | None:
    vals = []
    for i in (1,2):
        chrom = clean(r.get(f'CHROM_{i}')); pos = clean(r.get(f'POS_{i}'))
        ref = clean(r.get(f'REF_{i}')).upper(); alt = clean(r.get(f'ALT_{i}')).upper()
        if not all((chrom,pos,ref,alt)): return None
        if not chrom.lower().startswith('chr'): chrom = 'chr' + chrom
        vals.append((chrom, str(integer(pos,-1)), ref, alt))
    return tuple(sorted(vals))


def variant_display(r: dict[str, Any], i: int) -> str:
    return f"{clean(r.get(f'CHROM_{i}'))}:{clean(r.get(f'POS_{i}'))}:{clean(r.get(f'REF_{i}'))}>{clean(r.get(f'ALT_{i}'))}"


def explanation(r: dict[str, Any]) -> str:
    parts=[]
    if integer(r.get('PIVOT_ASSOCIATION_TIER')) == 0: parts.append('reviewed Mendelian gene-disease association')
    elif integer(r.get('PIVOT_ASSOCIATION_TIER')) == 2: parts.append('complex-susceptibility association; not Mendelian-pair support')
    else: parts.append('gene-disease association unclassified; review required')
    if integer(r.get('PIVOT_INHERITANCE_TIER')) == 0: parts.append('reviewed AR model with two heterozygous calls; phase unresolved')
    if integer(r.get('PIVOT_CLINVAR_DISEASE_MATCH_PLP')): parts.append('exact ClinVar P/LP evidence linked to target condition')
    elif integer(r.get('PIVOT_CLINVAR_PLP_EXACT')): parts.append('exact ClinVar P/LP exists, but disease-context match not established')
    if integer(r.get('PIVOT_HPO_DOMAIN_COUNT')): parts.append(f"phenotype support spans {integer(r.get('PIVOT_HPO_DOMAIN_COUNT'))} distinct domains")
    if integer(r.get('PIVOT_AF_RARE_LE_1E4')): parts.append(f"{integer(r.get('PIVOT_AF_RARE_LE_1E4'))} alleles have exact AF <=1e-4")
    if r.get('PIVOT_NEARBY_INDEL_REVIEW') == 'YES': parts.append('nearby indel(s): review representation/read support')
    if r.get('PIVOT_MHC_REGION_REVIEW') == 'YES': parts.append('MHC region: mapping/allele representation review')
    if integer(r.get('PIVOT_CLINVAR_BENIGN_EXACT')): parts.append('exact benign ClinVar allele present')
    if integer(r.get('PIVOT_CLINVAR_CONFLICTING_EXACT')): parts.append('conflicting ClinVar classification present')
    return '; '.join(parts)


def choose_draft(rows: list[dict[str, Any]], max_rows: int, max_pairs_per_gene: int, max_mhc: int) -> list[dict[str, Any]]:
    chosen=[]; seen_gene=Counter(); seen_pairs=set(); mhc_count=0
    for r in rows:
        gene=norm_gene(r.get('GENE',''))
        key=row_pair_key(r)
        if not key or key in seen_pairs: continue
        if seen_gene[gene] >= max_pairs_per_gene: continue
        if r.get('PIVOT_ASSOCIATION_CLASS') == 'COMPLEX_SUSCEPTIBILITY':
            # do not automatically exclude the gene/variants; exclude this pair from
            # the main Mendelian pair shortlist because association itself is complex.
            continue
        if r.get('PIVOT_MHC_REGION_REVIEW') == 'YES' and mhc_count >= max_mhc: continue
        chosen.append(r); seen_gene[gene]+=1; seen_pairs.add(key)
        if r.get('PIVOT_MHC_REGION_REVIEW') == 'YES': mhc_count += 1
        if len(chosen) >= max_rows: break
    return chosen


def write_submission(path: Path, chosen: list[dict[str, Any]], epcrs: list[float], primary_only: bool = False) -> None:
    cols=['proband_id','chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2','epcr','finding_type','notes']
    subset=chosen[:1] if primary_only else chosen
    with path.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fieldnames=cols, lineterminator='\n')
        w.writeheader()
        for i,r in enumerate(subset):
            epcr=epcrs[i] if i < len(epcrs) else max(0.01, epcrs[-1] - 0.03*(i-len(epcrs)+1))
            notes=explanation(r)
            if r.get('PIVOT_PHASE_STATUS'):
                notes += '; phase status: ' + clean(r.get('PIVOT_PHASE_STATUS'))
            w.writerow({
                'proband_id':'PROBAND01',
                'chrom_1': clean(r.get('CHROM_1')) if clean(r.get('CHROM_1')).startswith('chr') else 'chr'+clean(r.get('CHROM_1')),
                'pos_1': integer(r.get('POS_1'),0), 'ref_1':clean(r.get('REF_1')).upper(),'alt_1':clean(r.get('ALT_1')).upper(),
                'chrom_2': clean(r.get('CHROM_2')) if clean(r.get('CHROM_2')).startswith('chr') else 'chr'+clean(r.get('CHROM_2')),
                'pos_2': integer(r.get('POS_2'),0), 'ref_2':clean(r.get('REF_2')).upper(),'alt_2':clean(r.get('ALT_2')).upper(),
                'epcr': f'{epcr:.4f}', 'finding_type':'primary' if i==0 else 'secondary', 'notes':notes[:500]
            })


def main() -> int:
    ap=argparse.ArgumentParser(description='PIVOT-MVA: lightweight pair-level evidence ranking with phenotype-domain collapse and sensitivity analysis.')
    ap.add_argument('--pairs', default=str(DEFAULT_PAIRS), help='VariantForge v0.5.3 pair_candidates TSV')
    ap.add_argument('--hpo-summary', default=str(DEFAULT_HPO), help='Gene-disease HPO summary TSV with MATCHED_PHENOTYPES')
    ap.add_argument('--inheritance', default=str(DEFAULT_INHERITANCE), help='Curated inheritance TSV; optional')
    ap.add_argument('--out', default=str(DEFAULT_OUT), help='Output folder')
    ap.add_argument('--association-overrides', default=None, help='Reviewed disease-level association TSV; default created under output')
    ap.add_argument('--max-rows', type=int, default=10)
    ap.add_argument('--max-pairs-per-gene', type=int, default=1)
    ap.add_argument('--max-mhc-pairs', type=int, default=1)
    ap.add_argument('--primary-epcr', type=float, default=0.85, help='Judgement-based EPCR estimate, not calibrated')
    ap.add_argument('--secondary-epcr-start', type=float, default=0.45, help='Start for exploratory secondary rows; not calibrated')
    args=ap.parse_args()
    if not 1 <= args.max_rows <= 10: ap.error('--max-rows must be 1..10')
    if args.max_pairs_per_gene < 1: ap.error('--max-pairs-per-gene must be >=1')
    if args.max_mhc_pairs < 0: ap.error('--max-mhc-pairs must be >=0')
    if not 0 < args.primary_epcr <= 1: ap.error('--primary-epcr must be in (0,1]')
    if not 0 < args.secondary_epcr_start < args.primary_epcr: ap.error('--secondary-epcr-start must be in (0, primary-epcr)')
    pair_path=Path(args.pairs).expanduser(); hpo_path=Path(args.hpo_summary).expanduser()
    inh_path=Path(args.inheritance).expanduser() if args.inheritance else None
    out=Path(args.out).expanduser(); out.mkdir(parents=True, exist_ok=True)
    if not pair_path.exists(): raise FileNotFoundError(f'Pair TSV not found: {pair_path}')
    if not hpo_path.exists(): raise FileNotFoundError(f'HPO summary not found: {hpo_path}')
    association_path=Path(args.association_overrides).expanduser() if args.association_overrides else out/'pivot_mva_association_overrides.tsv'
    associations=load_associations(association_path)
    pairs,_=read_tsv(pair_path); hpo_rows,_=read_tsv(hpo_path)
    hpomap=hpo_index(hpo_rows); inheritance=load_inheritance(inh_path)
    if not pairs: raise RuntimeError('Pair input has no rows')
    work=[]
    for r in pairs:
        f=build_features(r,hpomap,associations,inheritance)
        f['PIVOT_EXPLANATION']=explanation(f)
        work.append(f)
    work.sort(key=rank_tuple)
    for i,r in enumerate(work,1): r['PIVOT_FINAL_RANK']=i
    write_tsv(out/'pivot_mva_scored_pairs.tsv',work)

    # Counterfactual sensitivity: independently suppress one phenotype domain
    # then rerank all pairs. This measures ranking dependence, not causality.
    top=work[:50]
    base_ranks={clean(r.get('PAIR_ID')):i for i,r in enumerate(work,1)}
    sens=[]
    for domain in DOMAIN_ORDER:
        ranked=sorted(work,key=lambda x: rank_tuple(x,omit_domain=domain))
        for i,r in enumerate(ranked,1):
            pid=clean(r.get('PAIR_ID')) or f"{r.get('GENE')}|{r.get('VARIANT_1_KEY')}|{r.get('VARIANT_2_KEY')}"
            if pid not in base_ranks: base_ranks[pid]=i
        rankmap={clean(r.get('PAIR_ID')) or f"{r.get('GENE')}|{r.get('VARIANT_1_KEY')}|{r.get('VARIANT_2_KEY')}":i for i,r in enumerate(ranked,1)}
        for r in top:
            pid=clean(r.get('PAIR_ID')) or f"{r.get('GENE')}|{r.get('VARIANT_1_KEY')}|{r.get('VARIANT_2_KEY')}"
            sensrow={
                'BASELINE_RANK':integer(r.get('PIVOT_FINAL_RANK')),'PAIR_ID':pid,'GENE':r.get('GENE',''),
                'OMITTED_DOMAIN':domain,'RANK_WITH_DOMAIN_OMITTED':rankmap.get(pid,'NA'),
                'BASELINE_DOMAIN_SCORE':r.get('PIVOT_HPO_DOMAIN_SCORE','0'),
                'DOMAIN_SCORE_REMOVED':r.get('PIVOT_HPO_'+domain.upper(),'0'),
                'INTERPRETATION':'counterfactual ranking sensitivity; not causal probability'
            }
            sens.append(sensrow)
    write_tsv(out/'pivot_mva_counterfactual_sensitivity.tsv',sens)

    # Main submission shortlist: one pair per gene; at most one chr6 MHC pair.
    chosen=choose_draft(work,args.max_rows,args.max_pairs_per_gene,args.max_mhc_pairs)
    if not chosen:
        raise RuntimeError('No candidates survived draft diversity filters')
    # One-row draft is primary recommendation. EPCR is a manual, non-calibrated
    # point estimate constrained to rank the primary first. Diverse file values are
    # explicitly marked judgement-based. Evaluator uses rank/threshold ordering.
    write_submission(out/'pivot_mva_submission_primary_only.csv',chosen,[args.primary_epcr],primary_only=True)
    epcrs=[args.primary_epcr]
    for i in range(1,len(chosen)):
        # lower secondary estimates; do not claim population calibration
        epcrs.append(max(0.01, args.secondary_epcr_start - 0.04*(i-1)))
    write_submission(out/'pivot_mva_submission_diverse10.csv',chosen,epcrs,primary_only=False)

    report=[]
    report += ['PIVOT-MVA v1.0 — Pairwise Inheritance, Phenotype-domain and Traceable evidence', '='*82,
        f'Input pair rows: {len(pairs)}',f'HPO disease associations available: {len(hpo_rows)}',
        f'Output directory: {out}',f'Association overrides: {association_path}',
        'Model type: deterministic, lexicographic evidence ranking; no trained ML and no probability calibration.',
        'No raw VCF is opened or reprocessed. Python standard library only.',
        '', 'MODEL DESIGN',
        '1. Disease-association gate: reviewed Mendelian relation > unclassified > complex susceptibility.',
        '2. Inheritance/genotype compatibility is assessed before phenotype tie-breaks; phase is never inferred from 0/1 GTs.',
        '3. Exact ClinVar condition-linked P/LP support is separate from unrelated P/LP, benign, conflicting and VUS evidence.',
        '4. HPO matches are grouped into five domains; within each domain only the most informative match contributes, limiting correlated-term double counting.',
        '5. Population AF is used only when present; missing AF is never converted to zero or rarity.',
        '6. DP/GQ, nearby indels and MHC location generate review flags; proximity or location do not prove artifacts.',
        '7. Counterfactual sensitivity reranks candidates after omitting each HPO domain.',
        '', 'PHENOTYPE DOMAIN IC VALUES (from existing local HPO report; relative support only)',
    ]
    for term,(domain,ic) in PHENOTYPES.items(): report.append(f'- {term}: domain={domain}; IC={ic}')
    report += ['', 'TOP 25 PAIRS (ordinal rank is not causal probability)']
    for r in work[:25]:
        report.append(
            f"{r['PIVOT_FINAL_RANK']}. {r.get('GENE','')} | disease={r.get('HPO_BEST_DISEASE_ID','')} | "
            f"assoc={r.get('PIVOT_ASSOCIATION_CLASS')} | inheritance={r.get('INHERITANCE_MODEL_HINT','UNKNOWN')} "
            f"({r.get('PIVOT_INHERITANCE_STATUS')}) | exact P/LP={r.get('PIVOT_CLINVAR_PLP_EXACT')} "
            f"(condition-linked={r.get('PIVOT_CLINVAR_DISEASE_MATCH_PLP')}) | HPO domains={r.get('PIVOT_HPO_DOMAIN_COUNT')} "
            f"(score={r.get('PIVOT_HPO_DOMAIN_SCORE')}) | rare AF count={r.get('PIVOT_AF_RARE_LE_1E4')} | "
            f"{variant_display(r,1)} + {variant_display(r,2)}")
        report.append('    rationale: '+clean(r.get('PIVOT_EXPLANATION','')))
    report += ['', 'SUBMISSION DRAFTS',
        f"Primary-only rows: 1; file=pivot_mva_submission_primary_only.csv; EPCR={args.primary_epcr:.4f} is a judgment-based, non-calibrated estimate.",
        f"Diverse exploratory rows: {len(chosen)}; file=pivot_mva_submission_diverse10.csv; secondary values are ordinal judgment estimates, not calibrated.",
        'The primary-only file is the precision-first option. The diverse file is for manual review, not automatic submission.',
        '', 'LIMITATIONS / REVIEW GATES',
        '- Exact ClinVar evidence is allele-specific; condition linkage is string-based for labels and must be manually reviewed.',
        '- The HPO evidence table is gene-disease annotation coverage, not a patient-specific likelihood model. Missing HPO matches are not contradictions.',
        '- The disease-association config is a small reviewed knowledge layer, not an exhaustive disease mechanism database.',
        '- No BAM/CRAM, trio or read-backed phasing is available through these inputs; BUB1B phase remains unresolved.',
        '- No independent calibration cohort is present. Never describe ordinal ranks or EPCR draft values as trained/calibrated probabilities.',
        '- Do not commit raw VCF, individual-level genotype tables or genome-wide variant tables to a public repo.',
        '', 'PUBLIC METHOD REFERENCES / CONTEXT',
        '- Official Track 1 evaluator: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py',
        '- Exomiser: https://github.com/exomiser/Exomiser (conceptual reference; not executed here).',
        '- LIRICAL: https://github.com/TheJacksonLaboratory/LIRICAL (existing local tool; its disease posterior is not EPCR).',
        '- ClinGen BUB1B/MVA1 gene-disease validity source recorded in association overrides.',
        '- This project does not copy third-party source code; public reports from participants are external corroboration only, not model inputs.'
    ]
    (out/'RUN_REPORT.txt').write_text('\n'.join(report)+'\n',encoding='utf-8')

    print('PIVOT-MVA v1.0 completed')
    print(f'Pair hypotheses scored: {len(work)}')
    print(f'HPO associations read: {len(hpo_rows)}')
    print(f'Output: {out}')
    print('Files: pivot_mva_scored_pairs.tsv, pivot_mva_counterfactual_sensitivity.tsv,')
    print('       pivot_mva_submission_primary_only.csv, pivot_mva_submission_diverse10.csv, RUN_REPORT.txt')
    print('\nTOP 15')
    for r in work[:15]:
        print(f"{integer(r.get('PIVOT_FINAL_RANK')):>4} | {clean(r.get('GENE')):<12} | assoc={clean(r.get('PIVOT_ASSOCIATION_CLASS')):<24} | inh={clean(r.get('INHERITANCE_MODEL_HINT')):<8} | P/LP={integer(r.get('PIVOT_CLINVAR_PLP_EXACT'))} condition-linked={integer(r.get('PIVOT_CLINVAR_DISEASE_MATCH_PLP'))} | HPO domains={integer(r.get('PIVOT_HPO_DOMAIN_COUNT'))} | {variant_display(r,1)} + {variant_display(r,2)}")
    print('\nPRIMARY-ONLY CSV (inspect manually before submission)')
    print((out/'pivot_mva_submission_primary_only.csv').read_text(encoding='utf-8'))
    return 0

if __name__ == '__main__':
    try: sys.exit(main())
    except Exception as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(2)
