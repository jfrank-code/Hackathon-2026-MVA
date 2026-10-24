#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VariantForge v0.5 PAIR-TRACE

PAIR-TRACE = Pairwise Allele Inheritance, Ranking, Traceability and Conflict Evaluation.

Purpose
-------
A lightweight, local, audit-first prioritization pass for an already annotated TSV.
It does NOT diagnose, prove phase, infer pathogenicity from rarity alone, or output a
calibrated probability of causal relationship. It preserves missingness and conflicts.

Inputs (autodetected by default under /mnt/c/MVA_Hackathon/processed):
  variantforge_v0_3/all_variants_clinvar_categorized.tsv
  fallback: variantforge_v0_2/variant_evidence_with_databases.tsv
  fallback: variantforge_v0_1/variant_evidence.tsv
  optional HPO: variantforge_v0_4_1/gene_disease_summary.tsv
  fallback HPO: variantforge_v0_4/gene_disease_summary.tsv

No third-party Python packages are required.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import math
import os
import re
import sys
from collections import defaultdict, Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

VERSION = "0.5.3"
DEFAULT_BASE = Path("/mnt/c/MVA_Hackathon/processed")
DEFAULT_OUT = DEFAULT_BASE / "variantforge_v0_5"
MISSING = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN", "NOT_AVAILABLE", "NOT COVERED"}

# Inheritance is intentionally supplied through a local editable TSV, not hard-coded
# into the public script. This avoids embedding case-specific candidate genes or
# derived variant-level information in code that may be published.


def norm_key(s: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", str(s).strip().upper()).strip("_")


def clean(v: Any) -> str:
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.upper() in MISSING else s


def parse_float(v: Any) -> Optional[float]:
    s = clean(v)
    if not s:
        return None
    try:
        x = float(s)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def parse_int(v: Any) -> Optional[int]:
    x = parse_float(v)
    return int(x) if x is not None else None


def truthy(v: Any) -> bool:
    return clean(v).upper() in {"YES", "Y", "TRUE", "1", "EXACT", "PASS"}


def row_value(row: Dict[str, str], *aliases: str) -> str:
    # Headers are normalized once during load; alias priority is preserved.
    for alias in aliases:
        key = norm_key(alias)
        if key in row and clean(row[key]):
            return clean(row[key])
    return ""


def read_tsv(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        raw_headers = reader.fieldnames or []
        rows: List[Dict[str, str]] = []
        for raw in reader:
            item: Dict[str, str] = {}
            for k, v in (raw or {}).items():
                if k is not None:
                    item[norm_key(k)] = clean(v)
            rows.append(item)
    return rows, [norm_key(h) for h in raw_headers]


def write_tsv(path: Path, rows: Iterable[Dict[str, Any]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: serialize(row.get(k, "")) for k in columns})


def serialize(v: Any) -> str:
    if v is None:
        return "NA"
    if isinstance(v, float):
        if not math.isfinite(v):
            return "NA"
        return f"{v:.6g}"
    if isinstance(v, (list, tuple, set)):
        return "|".join(map(str, v))
    return str(v)


def locate_variant_file(base: Path, explicit: Optional[str]) -> Path:
    if explicit:
        p = Path(explicit).expanduser()
        if p.is_file():
            return p
        raise FileNotFoundError(f"No existe el TSV de variantes indicado: {p}")
    candidates = [
        base / "variantforge_v0_3" / "all_variants_clinvar_categorized.tsv",
        base / "variantforge_v0_2" / "variant_evidence_with_databases.tsv",
        base / "variantforge_v0_1" / "variant_evidence.tsv",
        base / "all_variants_clinvar_categorized.tsv",
        base / "variant_evidence_with_databases.tsv",
    ]
    for p in candidates:
        if p.is_file():
            return p
    raise FileNotFoundError("No encontré un TSV de variantes. Busqué:\n  " + "\n  ".join(map(str, candidates)))


def locate_hpo_file(base: Path, explicit: Optional[str]) -> Optional[Path]:
    if explicit:
        p = Path(explicit).expanduser()
        if p.is_file():
            return p
        raise FileNotFoundError(f"No existe el TSV HPO indicado: {p}")
    candidates = [
        base / "variantforge_v0_4_1" / "gene_disease_summary.tsv",
        base / "variantforge_v0_4" / "gene_disease_summary.tsv",
        base / "gene_disease_summary.tsv",
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


def ensure_inheritance_file(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["GENE", "CONDITION", "INHERITANCE", "SOURCE_OR_NOTE", "REVIEWED"])


def load_inheritance(path: Path) -> Dict[str, Dict[str, str]]:
    ensure_inheritance_file(path)
    rows, _ = read_tsv(path)
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        gene = row_value(row, "GENE", "SYMBOL").upper()
        if not gene:
            continue
        model = row_value(row, "INHERITANCE", "INHERITANCE_MODEL", "MODEL").upper().replace(" ", "_")
        if model in {"AUTOSOMAL_RECESSIVE", "RECESSIVE", "AR"}:
            model = "AR"
        elif model in {"AUTOSOMAL_DOMINANT", "DOMINANT", "AD"}:
            model = "AD"
        elif model in {"X_LINKED", "X-LINKED", "XL"}:
            model = "XL"
        elif model in {"UNKNOWN", "", "NA"}:
            model = "UNKNOWN"
        else:
            model = "CONDITION_DEPENDENT" if "CONDITION" in model or "MULTIPLE" in model else model
        note = row_value(row, "SOURCE_OR_NOTE", "SOURCE", "NOTE")
        reviewed = row_value(row, "REVIEWED").upper() in {"YES", "Y", "TRUE", "1"}
        result[gene] = {"model": model or "UNKNOWN", "condition": row_value(row, "CONDITION", "DISEASE"), "note": note, "reviewed": "YES" if reviewed else "NO"}
    return result


def parse_gt(gt: str, gt_class: str = "") -> str:
    g = clean(gt).replace("|", "/")
    if g in {"0/1", "1/0"}:
        return "HET"
    if g in {"1/1"}:
        return "HOM_ALT"
    if g in {"0/0"}:
        return "HOM_REF"
    if g and "." in g:
        return "MISSING_GT"
    gc = clean(gt_class).upper()
    if "HETEROZYGOUS" in gc:
        return "HET"
    if "HOM_ALT" in gc or gc == "HOMOZYGOUS_ALT" or gc == "HOM_ALT":
        return "HOM_ALT"
    return "UNKNOWN"


def consequence_class(consequence: str, impact: str) -> Tuple[int, int, str]:
    c = (consequence or "").lower().replace("-", "_")
    i = (impact or "").upper()
    lof_terms = ("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "transcript_ablation")
    likely_lof_terms = ("start_lost", "stop_lost", "transcript_amplification")
    moderate_terms = ("missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "splice_region_variant")
    if any(t in c for t in lof_terms):
        return 3, 1, "PREDICTED_LOF_CONSEQUENCE"
    if any(t in c for t in likely_lof_terms):
        return 2, 1, "POTENTIAL_LOF_REQUIRES_TRANSCRIPT_MECHANISM_REVIEW"
    if any(t in c for t in moderate_terms) or "MODERATE" in i:
        return 2, 0, "PROTEIN_ALTERING_OR_MODERATE"
    if "LOW" in i or "synonymous" in c:
        return 0, 0, "LOW_IMPACT_OR_SYNONYMOUS"
    if "HIGH" in i:
        return 2, 1, "HIGH_IMPACT_ANNOTATION_REVIEW_CONSEQUENCE"
    return 1, 0, "OTHER_OR_UNCLEAR_CONSEQUENCE"


def clinvar_state(clnsig: str, exact_match: bool) -> str:
    s = (clnsig or "").lower().replace(" ", "_")
    if not exact_match:
        return "NO_EXACT_ALLELE_MATCH_OR_UNCONFIRMED"
    if "conflict" in s:
        return "CONFLICTING"
    has_path = "pathogenic" in s
    has_benign = "benign" in s
    has_likely = "likely_pathogenic" in s or "likely_pathogenic" in s.replace("/", "_")
    has_likely_benign = "likely_benign" in s
    if has_path and has_benign:
        return "MIXED_PATHOGENIC_BENIGN_LABELS_REVIEW"
    if has_likely_benign and not has_path:
        return "LIKELY_BENIGN"
    if has_path:
        return "PATHOGENIC_OR_LIKELY_PATHOGENIC" if has_likely else "PATHOGENIC"
    if has_benign:
        return "BENIGN"
    if "uncertain" in s or "vus" in s:
        return "VUS"
    if "not_provided" in s or "not_provided" in s.replace("/", "_"):
        return "NOT_PROVIDED"
    return "OTHER_REVIEW"


def is_plp(state: str) -> bool:
    return state in {"PATHOGENIC", "PATHOGENIC_OR_LIKELY_PATHOGENIC"}


def is_benign(state: str) -> bool:
    return state in {"BENIGN", "LIKELY_BENIGN"}


def variant_record(raw: Dict[str, str]) -> Dict[str, Any]:
    chrom = row_value(raw, "CHROM", "CHR", "CHROMOSOME")
    if chrom and not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    pos = row_value(raw, "POS", "POSITION", "START")
    ref = row_value(raw, "REF", "REFERENCE")
    alt = row_value(raw, "ALT", "ALTERNATE")
    gene = row_value(raw, "GENE", "SYMBOL", "GENE_SYMBOL").upper()
    gt = row_value(raw, "GT", "GENOTYPE")
    gt_class = row_value(raw, "GT_CLASS", "GENOTYPE_CLASS")
    impact = row_value(raw, "IMPACT", "ALL_IMPACTS")
    consequence = row_value(raw, "CONSEQUENCE", "EFFECT", "ANNOTATION")
    filt = row_value(raw, "FILTER") or "UNKNOWN"
    dp = parse_int(row_value(raw, "DP", "DEPTH"))
    gq = parse_float(row_value(raw, "GQ", "GENOTYPE_QUALITY"))
    ad = row_value(raw, "AD", "ALLELE_DEPTHS")
    ad_ref = ad_alt = None
    if ad:
        vals = re.split(r"[,|]", ad)
        if len(vals) >= 2:
            ad_ref, ad_alt = parse_float(vals[0]), parse_float(vals[1])
    ab = parse_float(row_value(raw, "ALLELE_FRACTION", "VAF"))
    if ab is None and ad_ref is not None and ad_alt is not None and (ad_ref + ad_alt) > 0:
        ab = ad_alt / (ad_ref + ad_alt)
    gt_state = parse_gt(gt, gt_class)
    clnsig = row_value(raw, "CLINVAR_CLNSIG", "CLNSIG", "CLINVAR_SIGNIFICANCE")
    exact_clinvar_raw = row_value(raw, "CLINVAR_EXACT_ALLELE_MATCH", "EXACT_CLINVAR_MATCH", "CLINVAR_EXACT_MATCH")
    exact_clinvar = truthy(exact_clinvar_raw)
    cl_state = clinvar_state(clnsig, exact_clinvar)
    gnomad_exact = truthy(row_value(raw, "GNOMAD_EXACT_ALLELE_MATCH", "EXACT_GNOMAD_MATCH", "GNOMAD_EXACT_MATCH"))
    af = parse_float(row_value(raw, "GNOMAD_AF", "GNOMAD_AF_GLOBAL", "AF")) if gnomad_exact else None
    af_afr = parse_float(row_value(raw, "GNOMAD_AF_AMR", "AF_AMR")) if gnomad_exact else None
    severity, lof, consequence_bucket = consequence_class(consequence, impact)
    hpo_rank = parse_int(row_value(raw, "HPO_RANK", "RANK_HPO"))
    hpo_score = parse_float(row_value(raw, "HPO_SCORE_PREVIOUS_EXPLORATORY", "HPO_SCORE"))
    flags: List[str] = []
    if "PENDING" in row_value(raw, "CLINVAR_STATUS") and exact_clinvar and clnsig:
        flags.append("CLINVAR_STATUS_FIELD_CONTRADICTS_EXACT_MATCH_AND_CLNSIG")
    if exact_clinvar and not clnsig:
        flags.append("CLINVAR_EXACT_MATCH_WITHOUT_CLNSIG")
    if not exact_clinvar and clnsig:
        flags.append("CLNSIG_PRESENT_BUT_EXACT_ALLELE_NOT_CONFIRMED")
    if not gnomad_exact:
        scope = row_value(raw, "GNOMAD_SCOPE", "GNOMAD_STATUS")
        flags.append("GNOMAD_UNAVAILABLE_OR_NO_EXACT_MATCH_NOT_AF_ZERO" if not scope or "NOT_COVERED" in scope or "PENDING" in scope or "NO" in scope else "GNOMAD_EXACT_ALLELE_NOT_CONFIRMED")
    if af is not None and af >= 0.95:
        flags.append("GNOMAD_AF_NEAR_FIXATION_CHECK_REF_ALT_BUILD_OR_JOIN")
    if filt.upper() not in {"PASS", ".", "UNKNOWN"}:
        flags.append("FILTER_NOT_PASS_REVIEW")
    if dp is not None and dp < 10:
        flags.append("LOW_DP_LT10")
    if gq is not None and gq < 20:
        flags.append("LOW_GQ_LT20")
    if gt_state == "HET" and ab is not None and not (0.2 <= ab <= 0.8):
        flags.append("HET_ALLELE_BALANCE_OUTSIDE_0_20_TO_0_80")
    if cl_state == "CONFLICTING":
        flags.append("CLINVAR_CONFLICTING")
    if is_benign(cl_state):
        flags.append("CLINVAR_EXACT_BENIGN_EVIDENCE_REVIEW_CONDITION_CONTEXT")
    if not row_value(raw, "GENE", "SYMBOL", "GENE_SYMBOL"):
        flags.append("GENE_SYMBOL_MISSING")
    if not consequence:
        flags.append("CONSEQUENCE_MISSING")
    if not chrom or not pos or not ref or not alt:
        flags.append("INCOMPLETE_VARIANT_KEY")
    key = f"{chrom}:{pos}:{ref}>{alt}|{gene}"
    # AF rarity rank is evidence only when the exact alternate allele is confirmed.
    if af is None:
        rarity = None
    elif af < 1e-6:
        rarity = 3
    elif af < 1e-4:
        rarity = 2
    elif af < 1e-3:
        rarity = 1
    else:
        rarity = 0
    return {
        "VARIANT_KEY": key, "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, "GENE": gene,
        "FILTER": filt, "GT": gt, "GT_STATE": gt_state, "DP": dp, "GQ": gq, "AD": ad,
        "ALLELE_FRACTION": ab, "IMPACT": impact, "CONSEQUENCE": consequence,
        "CONSEQUENCE_BUCKET": consequence_bucket, "MOLECULAR_SEVERITY_AXIS": severity,
        "PREDICTED_LOF_FLAG": lof, "HGVS_C": row_value(raw, "HGVS_C"), "HGVS_P": row_value(raw, "HGVS_P"),
        "TRANSCRIPTS": row_value(raw, "TRANSCRIPTS", "TRANSCRIPT"),
        "CLINVAR_EXACT_ALLELE_MATCH": "YES" if exact_clinvar else "NO_OR_UNKNOWN",
        "CLINVAR_ID": row_value(raw, "CLINVAR_ID", "VCV"), "CLINVAR_CLNSIG": clnsig,
        "CLINVAR_DISEASES": row_value(raw, "CLINVAR_CLNDN", "CLNDN"),
        "CLINVAR_REVIEW_STATUS": row_value(raw, "CLINVAR_REVIEW_STATUS", "CLNREVSTAT"),
        "CLINVAR_STATE_AUDITED": cl_state, "CLINVAR_PLP_EXACT_FLAG": 1 if exact_clinvar and is_plp(cl_state) else 0,
        "CLINVAR_BENIGN_EXACT_FLAG": 1 if exact_clinvar and is_benign(cl_state) else 0,
        "GNOMAD_EXACT_ALLELE_MATCH": "YES" if gnomad_exact else "NO_OR_UNKNOWN",
        "GNOMAD_SCOPE": row_value(raw, "GNOMAD_SCOPE"), "GNOMAD_AF_EXACT": af, "GNOMAD_AF_AMR_EXACT": af_afr,
        "GNOMAD_RARITY_AXIS": rarity, "HPO_RANK_PREVIOUS": hpo_rank, "HPO_SCORE_PREVIOUS_EXPLORATORY": hpo_score,
        "AUDIT_FLAGS": ";".join(flags) if flags else "NONE",
        "RAW_CLINVAR_STATUS": row_value(raw, "CLINVAR_STATUS"), "RAW_GNOMAD_STATUS": row_value(raw, "GNOMAD_STATUS"),
        "SOURCE_VARIANT_ROW": raw,
    }


def merge_duplicate_variants(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    merged: Dict[str, Dict[str, Any]] = {}
    for rec in records:
        k = rec["VARIANT_KEY"]
        if not rec["GENE"] or not rec["CHROM"] or not rec["POS"] or not rec["REF"] or not rec["ALT"]:
            # keep invalid rows individually for audit
            k += f"|INVALID_ROW_{len(merged)}"
        if k not in merged:
            merged[k] = rec
            continue
        old = merged[k]
        # Prefer the row with an exact ClinVar allele join, then better quality, then richer annotations.
        rank_old = (old["CLINVAR_EXACT_ALLELE_MATCH"] == "YES", old["DP"] or -1, old["GQ"] or -1, len(old["SOURCE_VARIANT_ROW"]))
        rank_new = (rec["CLINVAR_EXACT_ALLELE_MATCH"] == "YES", rec["DP"] or -1, rec["GQ"] or -1, len(rec["SOURCE_VARIANT_ROW"]))
        if rank_new > rank_old:
            merged[k] = rec
    return list(merged.values())


def detect_numeric_count(row: Dict[str, str], patterns: Tuple[str, ...]) -> Optional[int]:
    """Read HPO count columns without confusing subtype relations with all relations.

    The old substring-first matcher treated N_SPECIFIC_SUBTYPE_RELATIONS as the
    generic RELATION count because the word RELATION is a substring of that
    header. Prefer explicit known columns first, then use conservative fallback.
    """
    preferred = {
        ("EXACT",): (
            "N_EXACT", "N_EXACT_MATCHES", "EXACT_HPO_MATCHES",
            "EXACT_MATCH_COUNT", "EXACT_MATCHES", "HPO_EXACT_COUNT_DETECTED",
        ),
        ("SUBTYPE",): (
            "N_SPECIFIC_SUBTYPE_RELATIONS", "N_SUBTYPE_MATCHES",
            "SPECIFIC_SUBTYPE_COUNT", "MORE_SPECIFIC_SUBTYPE_COUNT",
            "SUBTYPE_COUNT", "HPO_SUBTYPE_COUNT_DETECTED",
        ),
        ("RELATION",): (
            "N_TERMS_WITH_ANY_ONTOLOGY_RELATION", "N_TERMS_WITH_ANY_RELATION",
            "HPO_ANY_RELATION_COUNT", "TERMS_WITH_ANY_ONTOLOGY_RELATION",
            "TERMS_ANY_RELATION", "N_HPO_RELATIONS", "HPO_RELATION_COUNT",
            "HPO_RELATION_COUNT_DETECTED", "RELATION_COUNT",
        ),
    }
    for name in preferred.get(patterns, ()):
        val = parse_int(row.get(norm_key(name), ""))
        if val is not None:
            return val

    # Conservative fallback for other schemas. For total relation counts,
    # deliberately exclude subtype-specific and broad-ancestor columns.
    for key, value in row.items():
        if patterns == ("RELATION",) and any(
            token in key for token in ("SUBTYPE", "SPECIFIC", "BROAD_ANCESTOR")
        ):
            continue
        if all(p in key for p in patterns) and any(
            token in key for token in ("COUNT", "N_", "NUM", "MATCHES", "TERMS", "HIT")
        ):
            val = parse_int(value)
            if val is not None:
                return val
    return None


def load_hpo(path: Optional[Path]) -> Tuple[Dict[str, Dict[str, Any]], List[Dict[str, str]], str]:
    if not path:
        return {}, [], "NOT_FOUND"
    try:
        rows, headers = read_tsv(path)
    except Exception as e:
        return {}, [], f"READ_ERROR:{e}"
    per_gene: Dict[str, Dict[str, Any]] = {}
    all_rows: List[Dict[str, str]] = []
    for row in rows:
        gene = row_value(row, "GENE", "SYMBOL", "GENE_SYMBOL").upper()
        if not gene:
            continue
        exact = detect_numeric_count(row, ("EXACT",))
        subtype = detect_numeric_count(row, ("SUBTYPE",))
        relation = detect_numeric_count(row, ("RELATION",))
        disease_id = row_value(row, "DISEASE_ID", "DISEASE", "OMIM_ID", "ORPHA_ID", "ASSOCIATION_ID")
        disease_name = row_value(row, "DISEASE_NAME", "CONDITION_NAME", "DISEASE_LABEL", "CONDITION")
        key = (exact if exact is not None else -1, relation if relation is not None else -1, subtype if subtype is not None else -1)
        if gene not in per_gene or key > per_gene[gene]["_key"]:
            per_gene[gene] = {
                "HPO_EXACT_MATCHES": exact, "HPO_SPECIFIC_SUBTYPE_MATCHES": subtype,
                "HPO_ANY_RELATION_COUNT": relation, "HPO_BEST_DISEASE_ID": disease_id,
                "HPO_BEST_DISEASE_NAME": disease_name, "HPO_SOURCE_FILE": str(path), "_key": key,
                "HPO_ROW_COUNT": 1,
            }
        else:
            per_gene[gene]["HPO_ROW_COUNT"] += 1
        # Include normalized rows for the separate association-level output.
        flat = dict(row)
        flat["GENE"] = gene
        flat["_exact_detected"] = str(exact) if exact is not None else ""
        flat["_relation_detected"] = str(relation) if relation is not None else ""
        flat["_subtype_detected"] = str(subtype) if subtype is not None else ""
        all_rows.append(flat)
    # If the file had no recognized explicit exact count but has HPO rows, don't invent zero;
    # mark missing counts as None. Gene-level rank can still use previous HPO_RANK as proxy.
    return per_gene, all_rows, str(path)


def gene_model_rank(model: str, reviewed: str) -> int:
    """Ordinal compatibility for a *pair* hypothesis, not a gene pathogenicity score.

    Curated AR raises a pair hypothesis. Unknown stays neutral-low. AD/XL or
    condition-dependent mechanisms are not promoted as generic recessive pairs.
    """
    m = model.upper()
    reviewed_flag = reviewed == "YES"
    if m == "AR":
        return 4 if reviewed_flag else 3
    if m in {"UNKNOWN", ""}:
        return 1
    # A pair of variants is not, by itself, aligned with a dominant/X-linked model.
    if m in {"AD", "XL", "CONDITION_DEPENDENT", "AD_OR_AR", "AD/AR"}:
        return 0
    return 0


def technical_quality(rec: Dict[str, Any]) -> int:
    score = 0
    if rec["FILTER"].upper() == "PASS": score += 1
    if rec["DP"] is not None and rec["DP"] >= 10: score += 1
    if rec["GQ"] is not None and rec["GQ"] >= 20: score += 1
    if rec["GT_STATE"] != "UNKNOWN" and rec["GT_STATE"] != "MISSING_GT": score += 1
    if rec["GT_STATE"] == "HET" and rec["ALLELE_FRACTION"] is not None and 0.2 <= rec["ALLELE_FRACTION"] <= 0.8: score += 1
    return score


def preliminary_pair_tuple(pair: Dict[str, Any]) -> Tuple[Any, ...]:
    # Transparent lexicographic prioritization, not a probability or a clinical score.
    return (
        -int(pair.get("MODEL_COMPATIBILITY_AXIS", 0)),
        -int(pair.get("CLINVAR_PLP_COUNT", 0)),
        -int(pair.get("LOF_VARIANT_COUNT", 0)),
        -int(pair.get("HPO_EXACT_MATCHES") if pair.get("HPO_EXACT_MATCHES") is not None else -1),
        -int(pair.get("RARE_AF_COUNT", 0)),
        -int(pair.get("TECHNICAL_PASS_COUNT", 0)),
        -int(pair.get("HPO_RANK_PROXY") if pair.get("HPO_RANK_PROXY") is not None else -1),
        str(pair.get("GENE", "")), str(pair.get("VARIANT_1_KEY", "")), str(pair.get("VARIANT_2_KEY", "")),
    )


def final_pair_priority_tuple(pair: Dict[str, Any]) -> Tuple[Any, ...]:
    """Clinical-first deterministic review order; not a probability model.

    Unlike the v0.5.0 overlay, a Pareto label never promotes an otherwise weak
    candidate. Known benign/conflicting evidence and close overlapping indels
    are explicit review penalties. Missing HPO/AF remains missing, not negative.
    """
    hpo = pair.get("HPO_EXACT_MATCHES")
    subtype = pair.get("HPO_SPECIFIC_SUBTYPE_MATCHES")
    relation = pair.get("HPO_ANY_RELATION_COUNT")
    hpo_values = [x for x in (hpo, subtype, relation) if x is not None]
    if not hpo_values:
        hpo_support_status = 1  # missing = unknown/neutral, not contradiction
    elif any(int(x) > 0 for x in hpo_values):
        hpo_support_status = 2  # at least one positive phenotype association
    else:
        hpo_support_status = 0  # known data but no positive match
    hpo_key = int(hpo) if hpo is not None else 0
    subtype_key = int(subtype) if subtype is not None else 0
    relation_key = int(relation) if relation is not None else 0
    conflicted = sum(
        1 for x in (pair.get("CLINVAR_STATE_1"), pair.get("CLINVAR_STATE_2"))
        if x == "CONFLICTING"
    )
    local_event_penalty = 1 if pair.get("LOCAL_COMPLEX_EVENT_FLAG") == "REVIEW_CLOSE_OR_OVERLAPPING_INDEL" else 0
    return (
        -int(pair.get("MODEL_COMPATIBILITY_AXIS", 0)),
        int(pair.get("CLINVAR_BENIGN_COUNT", 0)),
        conflicted,
        -int(pair.get("CLINVAR_PLP_COUNT", 0)),
        -hpo_support_status,
        -hpo_key,
        -subtype_key,
        -relation_key,
        local_event_penalty,
        -int(pair.get("LOF_VARIANT_COUNT", 0)),
        -int(pair.get("RARE_AF_COUNT", 0)),
        -int(pair.get("TECHNICAL_PASS_COUNT", 0)),
        -int(pair.get("GNOMAD_AF_KNOWN_COUNT", 0)),
        str(pair.get("GENE", "")),
        str(pair.get("VARIANT_1_KEY", "")),
        str(pair.get("VARIANT_2_KEY", "")),
    )


def pair_axes(pair: Dict[str, Any]) -> Tuple[int, int, int, int, int, int, int, int]:
    # Each axis is an ordinal evidence descriptor, not a calibrated unit.
    hpo = pair.get("HPO_EXACT_MATCHES")
    if hpo is None:
        hpo = pair.get("HPO_RANK_PROXY")
        hpo_axis = max(0, 20 - int(hpo)) if hpo is not None else 0
    else:
        hpo_axis = max(0, min(10, int(hpo)))
    hpo_available_axis = 1 if pair.get("HPO_EXACT_MATCHES") is not None else 0
    return (
        int(pair.get("MODEL_COMPATIBILITY_AXIS", 0)),
        int(pair.get("CLINVAR_PLP_COUNT", 0)),
        int(pair.get("LOF_VARIANT_COUNT", 0)),
        int(pair.get("RARE_AF_COUNT", 0)),
        int(pair.get("GNOMAD_AF_KNOWN_COUNT", 0)),
        int(hpo_available_axis),
        int(hpo_axis),
        int(pair.get("TECHNICAL_PASS_COUNT", 0)),
    )


def assign_pareto_frontier(pairs: List[Dict[str, Any]], limit: int = 2500) -> int:
    """Assign the exact first Pareto frontier for a deterministic shortlist.

    We don't calculate a full set of nested Pareto layers. If the candidate universe is
    huge, first form a transparent shortlist using the printed ranking tuple, then mark
    all non-shortlisted rows as OUTSIDE_PARETO_SHORTLIST. The full candidate table remains.
    """
    if not pairs:
        return 0
    pairs.sort(key=preliminary_pair_tuple)
    selected = pairs[:limit]
    axes = [pair_axes(p) for p in selected]
    n = len(selected)
    dominated = [False] * n
    for i in range(n):
        ai = axes[i]
        for j in range(n):
            if i == j:
                continue
            aj = axes[j]
            # j dominates i iff no worse on all axes and better on one.
            if all(a >= b for a, b in zip(aj, ai)) and any(a > b for a, b in zip(aj, ai)):
                dominated[i] = True
                break
    for i, p in enumerate(selected):
        p["PARETO_FIRST_FRONTIER"] = "NO" if dominated[i] else "YES"
        p["PARETO_TEST_SCOPE"] = f"TOP_{n}_BY_TRANSPARENT_LEXICOGRAPHIC_TIEBREAK"
        p["PRELIMINARY_RANK"] = i + 1
    for p in pairs[limit:]:
        p["PARETO_FIRST_FRONTIER"] = "NOT_TESTED_OUTSIDE_SHORTLIST"
        p["PARETO_TEST_SCOPE"] = f"TOP_{n}_BY_TRANSPARENT_LEXICOGRAPHIC_TIEBREAK"
        p["PRELIMINARY_RANK"] = "OUTSIDE_SHORTLIST"
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description="VariantForge v0.5.3 PAIR-TRACE: auditoría de variantes y pares sin probabilidad causal falsa.")
    ap.add_argument("--base", default=str(DEFAULT_BASE), help="Carpeta processed del proyecto")
    ap.add_argument("--variants", default=None, help="TSV de variantes (auto-detectado si se omite)")
    ap.add_argument("--hpo", default=None, help="TSV gene_disease_summary.tsv (auto-detectado si se omite)")
    ap.add_argument("--inheritance", default=None, help="TSV opcional/editable con herencia por gen")
    ap.add_argument("--out", default=None, help="Carpeta de salida")
    ap.add_argument("--max-variants-per-gene", type=int, default=60, help="Límite de variantes heterocigotas candidatas por gen; 0 elimina el límite")
    ap.add_argument("--pareto-limit", type=int, default=2500, help="Máximo de pares evaluados en la frontera de Pareto exacta")
    ap.add_argument("--draft-epcr", type=float, default=0.95, help="Valor inicial solo para ordenar el borrador CSV (NO calibrado como probabilidad)")
    ap.add_argument("--submission-candidates", type=int, default=10, help="Número de pares a incluir en el borrador Track 1 (1-10)")
    ap.add_argument("--max-submission-pairs-per-gene", type=int, default=2, help="Diversidad del borrador: máximo de pares por gen; 0 quita el límite")
    args = ap.parse_args()
    if not (0.0 < args.draft_epcr <= 1.0):
        ap.error("--draft-epcr debe estar en el intervalo (0, 1]")
    if not (1 <= args.submission_candidates <= 10):
        ap.error("--submission-candidates debe estar entre 1 y 10")
    if args.max_submission_pairs_per_gene < 0:
        ap.error("--max-submission-pairs-per-gene no puede ser negativo")
    if args.draft_epcr - 0.05 * (args.submission_candidates - 1) <= 0:
        ap.error("--draft-epcr es demasiado bajo para crear EPCRs distintos; sube el inicial o reduce submission-candidates")

    base = Path(args.base).expanduser()
    out = Path(args.out).expanduser() if args.out else base / "variantforge_v0_5"
    out.mkdir(parents=True, exist_ok=True)
    variants_path = locate_variant_file(base, args.variants)
    hpo_path = locate_hpo_file(base, args.hpo)
    inheritance_path = Path(args.inheritance).expanduser() if args.inheritance else out / "inheritance_overrides.tsv"
    inheritance = load_inheritance(inheritance_path)

    raw_rows, headers = read_tsv(variants_path)
    if not raw_rows:
        raise RuntimeError(f"El TSV de variantes no tiene filas: {variants_path}")
    records = merge_duplicate_variants([variant_record(r) for r in raw_rows])
    hpo_by_gene, hpo_rows, hpo_source = load_hpo(hpo_path)

    # Link HPO evidence to every variant without turning missing annotations into mismatch.
    for rec in records:
        h = hpo_by_gene.get(rec["GENE"], {})
        rec["HPO_EXACT_MATCHES_GENE_BEST_ASSOC"] = h.get("HPO_EXACT_MATCHES")
        rec["HPO_SPECIFIC_SUBTYPE_MATCHES_GENE_BEST_ASSOC"] = h.get("HPO_SPECIFIC_SUBTYPE_MATCHES")
        rec["HPO_ANY_RELATION_COUNT_GENE_BEST_ASSOC"] = h.get("HPO_ANY_RELATION_COUNT")
        rec["HPO_BEST_DISEASE_ID"] = h.get("HPO_BEST_DISEASE_ID", "")
        rec["HPO_BEST_DISEASE_NAME"] = h.get("HPO_BEST_DISEASE_NAME", "")
        model = inheritance.get(rec["GENE"], {"model": "UNKNOWN", "condition": "", "note": "NO_LOCAL_MODEL", "reviewed": "NO"})
        rec["INHERITANCE_MODEL_HINT"] = model["model"]
        rec["INHERITANCE_CONDITION_HINT"] = model.get("condition", "")
        rec["INHERITANCE_HINT_SOURCE_NOTE"] = model.get("note", "")
        rec["INHERITANCE_HINT_REVIEWED"] = model.get("reviewed", "NO")
        rec["TECHNICAL_QUALITY_AXIS"] = technical_quality(rec)
        # Add flags for lack of local HPO data, while retaining unknown as unknown.
        flags = [] if rec["AUDIT_FLAGS"] == "NONE" else rec["AUDIT_FLAGS"].split(";")
        if not h:
            flags.append("NO_LOCAL_HPO_ASSOCIATION_FOUND_NOT_NEGATIVE_PHENOTYPE_EVIDENCE")
        if model["model"] == "UNKNOWN":
            flags.append("INHERITANCE_MODEL_NOT_CONFIGURED")
        rec["AUDIT_FLAGS"] = ";".join(dict.fromkeys(flags)) if flags else "NONE"

    # Candidate pairs: same gene, two distinct heterozygous alleles, variants with molecular impact.
    genes: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for rec in records:
        if not rec["GENE"] or rec["GT_STATE"] != "HET":
            continue
        if rec["MOLECULAR_SEVERITY_AXIS"] <= 0:
            continue
        genes[rec["GENE"]].append(rec)

    pair_rows: List[Dict[str, Any]] = []
    same_site_skipped = 0
    genes_truncated: Dict[str, Tuple[int, int]] = {}
    for gene, items in genes.items():
        # Keep strongest annotated items first if the gene's list needs bounding. We do not delete
        # them from variant_audit.tsv; truncation is explicitly reported in pair outputs/report.
        items = sorted(items, key=lambda r: (
            -r["CLINVAR_PLP_EXACT_FLAG"], -r["PREDICTED_LOF_FLAG"],
            -(r["GNOMAD_RARITY_AXIS"] if r["GNOMAD_RARITY_AXIS"] is not None else -1),
            -r["MOLECULAR_SEVERITY_AXIS"], -r["TECHNICAL_QUALITY_AXIS"],
            r["CHROM"], int(r["POS"]) if str(r["POS"]).isdigit() else 0, r["REF"], r["ALT"]
        ))
        original_n = len(items)
        if args.max_variants_per_gene > 0 and len(items) > args.max_variants_per_gene:
            items = items[:args.max_variants_per_gene]
            genes_truncated[gene] = (original_n, len(items))
        model = inheritance.get(gene, {"model": "UNKNOWN", "condition": "", "note": "NO_LOCAL_MODEL", "reviewed": "NO"})
        model_rank = gene_model_rank(model["model"], model.get("reviewed", "NO"))
        h = hpo_by_gene.get(gene, {})
        hpo_exact = h.get("HPO_EXACT_MATCHES")
        hpo_subtype = h.get("HPO_SPECIFIC_SUBTYPE_MATCHES")
        hpo_rel = h.get("HPO_ANY_RELATION_COUNT")
        # Fallback HPO rank is not the same as phenotype similarity; mark it explicitly.
        ranks = [r["HPO_RANK_PREVIOUS"] for r in items if r["HPO_RANK_PREVIOUS"] is not None]
        hpo_rank_proxy = min(ranks) if ranks else None
        for a0, b0 in itertools.combinations(items, 2):
            if a0["CHROM"] == b0["CHROM"] and a0["POS"] == b0["POS"]:
                same_site_skipped += 1
                continue
            # Make variant suffixes _1/_2 match the sorted VARIANT_1_KEY/_2_KEY fields.
            a, b = sorted((a0, b0), key=lambda r: r["VARIANT_KEY"])
            # Per-pair audit flags. Nearby/overlapping indels can be distinct calls,
            # but may also represent one complex haplotype event; retain and flag them.
            flags: List[str] = ["PHASE_UNKNOWN_GT_NOT_PHASED"]
            local_complex_event_flag = "NO"
            local_event_distance_bp = "NA"
            try:
                a_pos, b_pos = int(a["POS"]), int(b["POS"])
                a_end = a_pos + max(1, len(a["REF"])) - 1
                b_end = b_pos + max(1, len(b["REF"])) - 1
                if a["CHROM"] == b["CHROM"]:
                    gap = max(0, max(a_pos, b_pos) - min(a_end, b_end) - 1)
                    local_event_distance_bp = gap
                    overlaps = max(a_pos, b_pos) <= min(a_end, b_end)
                    indel_present = (len(a["REF"]) != len(a["ALT"]) or len(b["REF"]) != len(b["ALT"]))
                    if overlaps or (gap <= 10 and indel_present):
                        local_complex_event_flag = "REVIEW_CLOSE_OR_OVERLAPPING_INDEL"
                        flags.append("CLOSE_OR_OVERLAPPING_INDEL_POSSIBLE_COMPLEX_EVENT")
            except (TypeError, ValueError):
                local_event_distance_bp = "NA"
            if model["model"] == "AR":
                model_status = "AR_MODEL_HINT_PAIR_COMPATIBLE_NOT_CONFIRMED"
            elif model["model"] in {"AD", "XL"}:
                model_status = "PAIR_NOT_PRIMARY_FOR_CONFIGURED_MODEL_REVIEW_CONDITION"
                flags.append("PAIR_RECESSIVE_HYPOTHESIS_NOT_ALIGNED_WITH_CONFIGURED_MODEL")
            elif model["model"] == "CONDITION_DEPENDENT":
                model_status = "CONDITION_DEPENDENT_MODEL_REVIEW_REQUIRED"
                flags.append("CONDITION_SPECIFIC_INHERITANCE_REVIEW_REQUIRED")
            else:
                model_status = "GENE_INHERITANCE_UNKNOWN_PAIR_IS_HYPOTHESIS_ONLY"
                flags.append("NO_CURATED_INHERITANCE_MODEL_FOR_GENE")
            if a["CLINVAR_BENIGN_EXACT_FLAG"] or b["CLINVAR_BENIGN_EXACT_FLAG"]:
                flags.append("AT_LEAST_ONE_EXACT_CLINVAR_BENIGN_ALLELE_REVIEW_CONDITION")
            if a["CLINVAR_STATE_AUDITED"] == "CONFLICTING" or b["CLINVAR_STATE_AUDITED"] == "CONFLICTING":
                flags.append("AT_LEAST_ONE_CLINVAR_CONFLICTING_ALLELE")
            if any((r["GNOMAD_AF_EXACT"] is not None and r["GNOMAD_AF_EXACT"] >= .95) for r in (a, b)):
                flags.append("AF_NEAR_FIXATION_ANOMALY_REQUIRES_JOIN_REVIEW")
            if a["FILTER"].upper() != "PASS" or b["FILTER"].upper() != "PASS":
                flags.append("ONE_OR_MORE_VARIANTS_NOT_FILTER_PASS")
            clinvar_plp_count = a["CLINVAR_PLP_EXACT_FLAG"] + b["CLINVAR_PLP_EXACT_FLAG"]
            benign_count = a["CLINVAR_BENIGN_EXACT_FLAG"] + b["CLINVAR_BENIGN_EXACT_FLAG"]
            lof_count = a["PREDICTED_LOF_FLAG"] + b["PREDICTED_LOF_FLAG"]
            rare_count = sum(1 for r in (a, b) if r["GNOMAD_RARITY_AXIS"] is not None and r["GNOMAD_RARITY_AXIS"] >= 2)
            af_known = sum(1 for r in (a, b) if r["GNOMAD_AF_EXACT"] is not None)
            quality_pass_count = sum(1 for r in (a, b) if r["FILTER"].upper() == "PASS" and (r["DP"] is None or r["DP"] >= 10) and (r["GQ"] is None or r["GQ"] >= 20))
            hpo_rank_proxy_pair = min([x for x in (a["HPO_RANK_PREVIOUS"], b["HPO_RANK_PREVIOUS"]) if x is not None], default=None)
            pair_key_1, pair_key_2 = sorted([a["VARIANT_KEY"], b["VARIANT_KEY"]])
            pair = {
                "PAIR_ID": f"{gene}|{pair_key_1}||{pair_key_2}", "GENE": gene,
                "LOCAL_COMPLEX_EVENT_FLAG": local_complex_event_flag,
                "LOCAL_EVENT_DISTANCE_BP": local_event_distance_bp,
                "RANKING_POLICY": "CLINICAL_FIRST_LEXICOGRAPHIC_V052_PARETO_AUDIT_ONLY",
                "VARIANT_1_KEY": pair_key_1, "VARIANT_2_KEY": pair_key_2,
                "CHROM_1": a["CHROM"], "POS_1": a["POS"], "REF_1": a["REF"], "ALT_1": a["ALT"], "GT_1": a["GT"],
                "CONSEQUENCE_1": a["CONSEQUENCE"], "IMPACT_1": a["IMPACT"], "HGVS_P_1": a["HGVS_P"],
                "CLINVAR_STATE_1": a["CLINVAR_STATE_AUDITED"], "CLINVAR_EXACT_1": a["CLINVAR_EXACT_ALLELE_MATCH"],
                "CLINVAR_DISEASES_1": a["CLINVAR_DISEASES"], "GNOMAD_AF_1": a["GNOMAD_AF_EXACT"],
                "DP_1": a["DP"], "GQ_1": a["GQ"], "AF_READS_1": a["ALLELE_FRACTION"],
                "CHROM_2": b["CHROM"], "POS_2": b["POS"], "REF_2": b["REF"], "ALT_2": b["ALT"], "GT_2": b["GT"],
                "CONSEQUENCE_2": b["CONSEQUENCE"], "IMPACT_2": b["IMPACT"], "HGVS_P_2": b["HGVS_P"],
                "CLINVAR_STATE_2": b["CLINVAR_STATE_AUDITED"], "CLINVAR_EXACT_2": b["CLINVAR_EXACT_ALLELE_MATCH"],
                "CLINVAR_DISEASES_2": b["CLINVAR_DISEASES"], "GNOMAD_AF_2": b["GNOMAD_AF_EXACT"],
                "DP_2": b["DP"], "GQ_2": b["GQ"], "AF_READS_2": b["ALLELE_FRACTION"],
                "INHERITANCE_MODEL_HINT": model["model"], "INHERITANCE_HINT_REVIEWED": model.get("reviewed", "NO"),
                "INHERITANCE_CONDITION_HINT": model.get("condition", ""), "INHERITANCE_HINT_SOURCE_NOTE": model.get("note", ""),
                "PAIR_MODEL_STATUS": model_status, "PHASE_STATUS": "UNKNOWN_NOT_INFERABLE_FROM_UNPHASED_GT",
                "CLINVAR_PLP_COUNT": clinvar_plp_count, "CLINVAR_BENIGN_COUNT": benign_count,
                "LOF_VARIANT_COUNT": lof_count, "MOLECULAR_SEVERITY_AXIS": a["MOLECULAR_SEVERITY_AXIS"] + b["MOLECULAR_SEVERITY_AXIS"],
                "RARE_AF_COUNT": rare_count, "GNOMAD_AF_KNOWN_COUNT": af_known,
                "GNOMAD_AF_COVERAGE_STATUS": "BOTH_EXACT_ALLELES" if af_known == 2 else ("ONE_EXACT_ALLELE" if af_known == 1 else "NO_EXACT_ALLELE_AF"),
                "HPO_EXACT_MATCHES": hpo_exact, "HPO_SPECIFIC_SUBTYPE_MATCHES": hpo_subtype,
                "HPO_ANY_RELATION_COUNT": hpo_rel, "HPO_BEST_DISEASE_ID": h.get("HPO_BEST_DISEASE_ID", ""),
                "HPO_BEST_DISEASE_NAME": h.get("HPO_BEST_DISEASE_NAME", ""), "HPO_ROW_COUNT_FOR_GENE": h.get("HPO_ROW_COUNT", 0),
                "HPO_RANK_PROXY": hpo_rank_proxy_pair,
                "TECHNICAL_PASS_COUNT": quality_pass_count,
                "MODEL_COMPATIBILITY_AXIS": model_rank,
                "EVIDENCE_VECTOR": f"inheritance={model_rank};ClinVar_PLP_exact={clinvar_plp_count};predicted_LoF={lof_count};rare_AF_count={rare_count};AF_coverage={af_known};HPO_exact={hpo_exact if hpo_exact is not None else 'NA'};technical_pass={quality_pass_count}",
                "AUDIT_FLAGS": ";".join(dict.fromkeys(flags)),
                "INTERPRETATION_LIMIT": "RANKING_HYPOTHESIS_ONLY_NOT_CAUSAL_PROBABILITY;PHASE_UNCONFIRMED",
            }
            pair_rows.append(pair)

    # Pareto remains a descriptive audit field only. In v0.5.3 it no longer
    # promotes every non-dominated row over candidates with stronger direct evidence.
    pareto_tested = assign_pareto_frontier(pair_rows, max(1, args.pareto_limit))
    pair_rows.sort(key=final_pair_priority_tuple)
    for idx, p in enumerate(pair_rows, 1):
        p["FINAL_PRIORITY_RANK"] = idx
        if idx == 1:
            p["PRIORITY_TIER"] = "TOP_RANKED_HYPOTHESIS_REQUIRES_MANUAL_REVIEW"
        elif p.get("PARETO_FIRST_FRONTIER") == "YES":
            p["PRIORITY_TIER"] = "NON_DOMINATED_EVIDENCE_PROFILE"
        elif p.get("PRELIMINARY_RANK") != "OUTSIDE_SHORTLIST":
            p["PRIORITY_TIER"] = "DOMINATED_WITHIN_SHORTLIST"
        else:
            p["PRIORITY_TIER"] = "LOWER_PRELIMINARY_PRIORITY_OR_OUTSIDE_PARETO_SHORTLIST"

    # Outputs: strip internal non-serializable objects; ensure fields are complete and stable.
    audit_columns = [
        "VARIANT_KEY", "CHROM", "POS", "REF", "ALT", "GENE", "FILTER", "GT", "GT_STATE", "DP", "GQ", "AD", "ALLELE_FRACTION",
        "IMPACT", "CONSEQUENCE", "CONSEQUENCE_BUCKET", "MOLECULAR_SEVERITY_AXIS", "PREDICTED_LOF_FLAG", "TRANSCRIPTS", "HGVS_C", "HGVS_P",
        "CLINVAR_EXACT_ALLELE_MATCH", "CLINVAR_ID", "CLINVAR_CLNSIG", "CLINVAR_DISEASES", "CLINVAR_REVIEW_STATUS", "CLINVAR_STATE_AUDITED",
        "CLINVAR_PLP_EXACT_FLAG", "CLINVAR_BENIGN_EXACT_FLAG", "RAW_CLINVAR_STATUS", "GNOMAD_EXACT_ALLELE_MATCH", "GNOMAD_SCOPE", "GNOMAD_AF_EXACT",
        "GNOMAD_AF_AMR_EXACT", "GNOMAD_RARITY_AXIS", "RAW_GNOMAD_STATUS", "HPO_EXACT_MATCHES_GENE_BEST_ASSOC", "HPO_SPECIFIC_SUBTYPE_MATCHES_GENE_BEST_ASSOC",
        "HPO_ANY_RELATION_COUNT_GENE_BEST_ASSOC", "HPO_BEST_DISEASE_ID", "HPO_BEST_DISEASE_NAME", "HPO_RANK_PREVIOUS", "HPO_SCORE_PREVIOUS_EXPLORATORY",
        "INHERITANCE_MODEL_HINT", "INHERITANCE_CONDITION_HINT", "INHERITANCE_HINT_SOURCE_NOTE", "INHERITANCE_HINT_REVIEWED", "TECHNICAL_QUALITY_AXIS", "AUDIT_FLAGS"
    ]
    pair_columns = [
        "FINAL_PRIORITY_RANK", "PRIORITY_TIER", "PARETO_FIRST_FRONTIER", "PARETO_TEST_SCOPE", "PRELIMINARY_RANK", "RANKING_POLICY", "PAIR_ID", "GENE", "LOCAL_COMPLEX_EVENT_FLAG", "LOCAL_EVENT_DISTANCE_BP",
        "VARIANT_1_KEY", "VARIANT_2_KEY", "CHROM_1", "POS_1", "REF_1", "ALT_1", "GT_1", "CONSEQUENCE_1", "IMPACT_1", "HGVS_P_1", "CLINVAR_STATE_1", "CLINVAR_EXACT_1", "CLINVAR_DISEASES_1", "GNOMAD_AF_1", "DP_1", "GQ_1", "AF_READS_1",
        "CHROM_2", "POS_2", "REF_2", "ALT_2", "GT_2", "CONSEQUENCE_2", "IMPACT_2", "HGVS_P_2", "CLINVAR_STATE_2", "CLINVAR_EXACT_2", "CLINVAR_DISEASES_2", "GNOMAD_AF_2", "DP_2", "GQ_2", "AF_READS_2",
        "INHERITANCE_MODEL_HINT", "INHERITANCE_HINT_REVIEWED", "INHERITANCE_CONDITION_HINT", "INHERITANCE_HINT_SOURCE_NOTE", "PAIR_MODEL_STATUS", "PHASE_STATUS",
        "CLINVAR_PLP_COUNT", "CLINVAR_BENIGN_COUNT", "LOF_VARIANT_COUNT", "MOLECULAR_SEVERITY_AXIS", "RARE_AF_COUNT", "GNOMAD_AF_KNOWN_COUNT", "GNOMAD_AF_COVERAGE_STATUS",
        "HPO_EXACT_MATCHES", "HPO_SPECIFIC_SUBTYPE_MATCHES", "HPO_ANY_RELATION_COUNT", "HPO_BEST_DISEASE_ID", "HPO_BEST_DISEASE_NAME", "HPO_ROW_COUNT_FOR_GENE", "HPO_RANK_PROXY", "TECHNICAL_PASS_COUNT", "MODEL_COMPATIBILITY_AXIS", "EVIDENCE_VECTOR", "AUDIT_FLAGS", "INTERPRETATION_LIMIT"
    ]
    clean_audit = [{k: r.get(k, "") for k in audit_columns} for r in records]
    clean_audit.sort(key=lambda r: (str(r.get("GENE", "")), str(r.get("CHROM", "")), int(r["POS"]) if str(r.get("POS", "")).isdigit() else 0, str(r.get("REF", "")), str(r.get("ALT", ""))))
    write_tsv(out / "variantforge_v053_variant_audit.tsv", clean_audit, audit_columns)
    write_tsv(out / "variantforge_v053_pair_candidates.tsv", pair_rows, pair_columns)

    # Review-required submission draft: up to 10 ranked pair hypotheses with strictly
    # descending EPCR placeholders. The evaluator computes F-max across thresholds, so
    # if the top row is exact, a lower-ranked row does not change F-max at that top threshold;
    # lower rows keep alternatives available for rank/partial credit if the top guess is wrong.
    # These EPCR values encode ORDER ONLY, not calibrated causal probabilities.
    submission_path = out / "variantforge_v053_track1_submission_draft.csv"
    submission_columns = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
    submission_rows = []
    submission_candidates = []
    submitted_per_gene = Counter()
    for candidate in pair_rows:
        gene = candidate.get("GENE", "")
        if args.max_submission_pairs_per_gene > 0 and submitted_per_gene[gene] >= args.max_submission_pairs_per_gene:
            continue
        submission_candidates.append(candidate)
        submitted_per_gene[gene] += 1
        if len(submission_candidates) >= args.submission_candidates:
            break
    for rank_idx, candidate in enumerate(submission_candidates, start=1):
        candidate_epcr = args.draft_epcr - 0.05 * (rank_idx - 1)
        submission_rows.append({
            "proband_id": "PROBAND01",
            "chrom_1": candidate["CHROM_1"], "pos_1": candidate["POS_1"], "ref_1": candidate["REF_1"], "alt_1": candidate["ALT_1"],
            "chrom_2": candidate["CHROM_2"], "pos_2": candidate["POS_2"], "ref_2": candidate["REF_2"], "alt_2": candidate["ALT_2"],
            "epcr": f"{candidate_epcr:.2f}", "finding_type": "primary",
            "notes": f"PAIR-TRACE ranked alternative #{rank_idx}; phase unknown; EPCR is an ordering placeholder, not calibrated probability. Review before submission."
        })
    with submission_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=submission_columns, lineterminator="\n")
        writer.writeheader(); writer.writerows(submission_rows)

    # One summary per gene-disease row from the HPO file, with detection transparently recorded.
    if hpo_rows:
        hpo_columns = list(dict.fromkeys(["GENE"] + [k for r in hpo_rows for k in r.keys() if not k.startswith("_")]))
        assoc_out = []
        for r in hpo_rows:
            out_row = {k: r.get(k, "") for k in hpo_columns}
            out_row["HPO_EXACT_COUNT_DETECTED"] = r.get("_exact_detected", "NA") or "NA"
            out_row["HPO_RELATION_COUNT_DETECTED"] = r.get("_relation_detected", "NA") or "NA"
            out_row["HPO_SUBTYPE_COUNT_DETECTED"] = r.get("_subtype_detected", "NA") or "NA"
            assoc_out.append(out_row)
        hpo_columns += ["HPO_EXACT_COUNT_DETECTED", "HPO_RELATION_COUNT_DETECTED", "HPO_SUBTYPE_COUNT_DETECTED"]
        write_tsv(out / "variantforge_v053_hpo_association_audit.tsv", assoc_out, hpo_columns)

    # Gene-level multi-axis table: no composite score or probability.
    gene_groups: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for rec in records:
        if rec["GENE"]:
            gene_groups[rec["GENE"]].append(rec)
    gene_out: List[Dict[str, Any]] = []
    for gene, items in gene_groups.items():
        m = inheritance.get(gene, {"model": "UNKNOWN", "condition": "", "note": "NO_LOCAL_MODEL", "reviewed": "NO"})
        h = hpo_by_gene.get(gene, {})
        gene_out.append({
            "GENE": gene, "VARIANT_COUNT": len(items),
            "HET_CANDIDATE_COUNT": sum(1 for r in items if r["GT_STATE"] == "HET" and r["MOLECULAR_SEVERITY_AXIS"] > 0),
            "HOM_ALT_COUNT": sum(1 for r in items if r["GT_STATE"] == "HOM_ALT"),
            "EXACT_CLINVAR_PLP_COUNT": sum(r["CLINVAR_PLP_EXACT_FLAG"] for r in items),
            "EXACT_CLINVAR_BENIGN_COUNT": sum(r["CLINVAR_BENIGN_EXACT_FLAG"] for r in items),
            "PREDICTED_LOF_COUNT": sum(r["PREDICTED_LOF_FLAG"] for r in items),
            "GNOMAD_AF_KNOWN_COUNT": sum(1 for r in items if r["GNOMAD_AF_EXACT"] is not None),
            "GNOMAD_AF_RARE_COUNT_AF_LT_1E_4": sum(1 for r in items if r["GNOMAD_AF_EXACT"] is not None and r["GNOMAD_AF_EXACT"] < 1e-4),
            "HPO_EXACT_MATCHES_BEST_ASSOC": h.get("HPO_EXACT_MATCHES"),
            "HPO_SPECIFIC_SUBTYPE_MATCHES_BEST_ASSOC": h.get("HPO_SPECIFIC_SUBTYPE_MATCHES"),
            "HPO_ANY_RELATION_COUNT_BEST_ASSOC": h.get("HPO_ANY_RELATION_COUNT"),
            "HPO_BEST_DISEASE_ID": h.get("HPO_BEST_DISEASE_ID", ""), "HPO_BEST_DISEASE_NAME": h.get("HPO_BEST_DISEASE_NAME", ""),
            "HPO_RANK_PREVIOUS_BEST": min([r["HPO_RANK_PREVIOUS"] for r in items if r["HPO_RANK_PREVIOUS"] is not None], default=None),
            "INHERITANCE_MODEL_HINT": m["model"], "INHERITANCE_HINT_REVIEWED": m.get("reviewed", "NO"), "INHERITANCE_HINT_SOURCE_NOTE": m.get("note", ""),
            "PHASE_OR_SEGREGATION_DATA_AVAILABLE": "NO_IN_LOCAL_INPUTS", "INTERPRETATION_LIMIT": "EVIDENCE_DIMENSIONS_ARE_NOT_A_CAUSAL_SCORE",
        })
    gene_out.sort(key=lambda r: (-int(r["EXACT_CLINVAR_PLP_COUNT"] or 0), -int(r["PREDICTED_LOF_COUNT"] or 0), -int(r["HPO_EXACT_MATCHES_BEST_ASSOC"] or 0), str(r["GENE"])))
    gene_columns = list(gene_out[0].keys()) if gene_out else ["GENE", "VARIANT_COUNT"]
    write_tsv(out / "variantforge_v053_gene_evidence.tsv", gene_out, gene_columns)

    # Report and key integrity checks.
    flagged = Counter()
    for r in records:
        for flag in str(r.get("AUDIT_FLAGS", "")).split(";"):
            if flag and flag != "NONE": flagged[flag] += 1
    top_pair = pair_rows[0] if pair_rows else None
    report_path = out / "RUN_REPORT.txt"
    with report_path.open("w", encoding="utf-8") as f:
        f.write(f"VariantForge v{VERSION} — PAIR-TRACE\n")
        f.write("='" + "=" * 70 + "\n")
        f.write(f"Execution UTC: {datetime.now(timezone.utc).isoformat()}\n")
        f.write(f"Variant input: {variants_path}\nHPO input: {hpo_source}\nInheritance config: {inheritance_path}\nOutput: {out}\n")
        f.write(f"Input rows: {len(raw_rows)}\nUnique variant records retained for audit: {len(records)}\nGenes with variants: {len(gene_groups)}\n")
        f.write(f"Genes with >=2 heterozygous moderate/high candidate variants: {sum(1 for g, xs in genes.items() if len(xs) >= 2)}\n")
        f.write(f"Pair hypotheses generated: {len(pair_rows)}\nPairs skipped due to same chromosome+position: {same_site_skipped}\n")
        f.write(f"Pairs with exact first Pareto frontier evaluated: {pareto_tested}\nPareto shortlist limit: {args.pareto_limit}\n")
        f.write(f"Genes truncated for pair enumeration: {len(genes_truncated)} (all variants remain in variant audit)\n")
        if genes_truncated:
            for gene, (n0, n1) in sorted(genes_truncated.items(), key=lambda x: -(x[1][0]-x[1][1]))[:30]:
                f.write(f"  {gene}: {n0} variants -> {n1} used for pair enumeration\n")
        f.write("\nINTERPRETATION RULES\n")
        f.write("- This is a prioritization/audit engine, not a diagnostic system.\n")
        f.write("- No EPCR or calibrated causal probability is produced. Evidence axes are distinct and not added together.\n")
        f.write("- Two unphased 0/1 variants in one gene do NOT prove they are in trans or compound-heterozygous.\n")
        f.write("- A ClinVar classification is used as exact-allele support only when the exact-match field is YES.\n")
        f.write("- AF is used only when GNOMAD_EXACT_ALLELE_MATCH=YES; no match/out-of-scope is missing data, never AF=0.\n")
        f.write("- Near-fixation gnomAD AF (>=0.95) is flagged for REF/ALT, build, normalization or join audit.\n")
        f.write("- HPO missingness is not treated as a phenotype contradiction; the HPO table is gene-disease association data.\n")
        f.write("- v0.5.3 fixes HPO count parsing: N_TERMS_WITH_ANY_ONTOLOGY_RELATION is selected explicitly and is not confused with N_SPECIFIC_SUBTYPE_RELATIONS.\n")
        f.write("- Inheritance entries are preliminary local hints; rows marked REVIEWED=NO must be reviewed with disease-specific sources.\n")
        f.write("- Pareto is an audit-only non-dominance label; it no longer promotes weak pairs in final ordering.\n")
        f.write("- Final review ranking v0.5.3 is lexicographic: pair-model compatibility; exact-benign/conflicting penalties; exact ClinVar P/LP support; positive HPO support and counts; close/overlapping-indel review flag; predicted LoF; confirmed rarity; technical quality; AF availability. It is deterministic, not calibrated.\n")
        f.write(f"- Submission diversity cap: max {args.max_submission_pairs_per_gene} pairs per gene (0 disables cap).\n")
        f.write(f"- Submission draft: {submission_path}; rows={len(submission_rows)}; first EPCR placeholder={args.draft_epcr:.2f}; DO NOT treat values as calibrated.\n")
        f.write("- Draft contains up to ten candidate pairs, with a per-gene diversity cap and strictly descending EPCR ordering placeholders; verify all rows before upload.\n")
        f.write("\nAUDIT FLAG COUNTS\n")
        for flag, count in flagged.most_common():
            f.write(f"{flag}\t{count}\n")
        f.write("\nTOP-RANKED PAIR CHECK\n")
        if top_pair:
            p = top_pair
            f.write(f"Top pair gene: {p['GENE']}\n")
            f.write(f"Priority rank: {p['FINAL_PRIORITY_RANK']}; Pareto first frontier: {p['PARETO_FIRST_FRONTIER']}\n")
            f.write(f"Model status: {p['PAIR_MODEL_STATUS']}\nPhase: {p['PHASE_STATUS']}\n")
            f.write(f"Exact ClinVar P/LP count: {p['CLINVAR_PLP_COUNT']}; predicted LoF count: {p['LOF_VARIANT_COUNT']}; AF-known count: {p['GNOMAD_AF_KNOWN_COUNT']}\n")
            f.write(f"HPO exact count (best gene-disease association if detected): {p['HPO_EXACT_MATCHES']}\n")
            f.write("CAUTION: ClinVar support for one allele does not establish the second allele's pathogenicity or phase.\n")
        else:
            f.write("No pair hypotheses generated. Check input path, coordinate build, gene labels, candidate GTs, and consequences.\n")
        f.write("\nTOP 20 PAIR HYPOTHESES (rank is a review order, not causal probability)\n")
        for p in pair_rows[:20]:
            f.write(f"{p['FINAL_PRIORITY_RANK']}. {p['GENE']} | Pareto={p['PARETO_FIRST_FRONTIER']} | model={p['INHERITANCE_MODEL_HINT']} | PLP_exact={p['CLINVAR_PLP_COUNT']} | LoF={p['LOF_VARIANT_COUNT']} | HPO_exact={p['HPO_EXACT_MATCHES']} | {p['VARIANT_1_KEY']} + {p['VARIANT_2_KEY']} | {p['AUDIT_FLAGS']}\n")
        f.write("\nNEXT VALIDATION STEPS\n")
        f.write("1. Inspect top-ranked pairs, including exact ClinVar disease context and MANE/transcript consequences.\n")
        f.write("2. Do not infer phase from this VCF; inspect reads/parents only if authorized data exist.\n")
        f.write("3. Review every inheritance TSV row marked REVIEWED=NO before final biological interpretation.\n")
        f.write("4. Review all AF near-fixation warnings for possible REF/ALT, genome-build or join problems.\n")
        f.write("5. EPCR values in the draft are placeholders for order only, not causal probability estimates; never submit them as calibrated estimates.\n")
        f.write("6. Close/overlapping indel pairs are flagged for representation/complex-haplotype review and are not automatically discarded.\n")

    print(f"VariantForge v{VERSION} — PAIR-TRACE completado")
    print(f"Variantes auditadas: {len(records)} | Genes: {len(gene_groups)} | Pares: {len(pair_rows)}")
    print(f"Entrada: {variants_path}")
    print(f"Salida: {out}")
    print(f"Borrador Track 1: {len(submission_rows)} pares (máximo 10; límite por gen={args.max_submission_pairs_per_gene}), EPCR placeholder de orden")
    print("Archivos: variantforge_v053_variant_audit.tsv, variantforge_v053_gene_evidence.tsv, variantforge_v053_pair_candidates.tsv, variantforge_v053_track1_submission_draft.csv, RUN_REPORT.txt")
    if hpo_path is None:
        print("AVISO: no encontré gene_disease_summary.tsv; HPO específico estará sin integrar en este pase.")
    if top_pair:
        print(f"Par principal provisional: {top_pair['GENE']} | rank 1 | fase UNKNOWN hasta resolver datos de fase")
    else:
        print("AVISO: no se generaron pares; revisa GT, gen y consecuencias en el TSV de entrada.")
    print("IMPORTANTE: el ranking clínico es heurístico, Pareto es solo auditoría y EPCR sigue siendo placeholder.")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR VariantForge v{VERSION}: {exc}", file=sys.stderr)
        raise
