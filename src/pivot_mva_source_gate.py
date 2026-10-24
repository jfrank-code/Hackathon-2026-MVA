#!/usr/bin/env python3
"""PIVOT-MVA Source Concordance Gate.

Locally checks every submitted allele against the original VCF using bcftools.
It does not infer pathogenicity. It audits exact CHROM/POS/REF/ALT identity and,
where available, the selected sample's GT/DP/GQ/AD. Output can contain sensitive
sample-level genotypes; keep it out of public repositories.

Requires Python 3.10+ and bcftools. No Python packages are needed.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

REQUIRED = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
            "chrom_2", "pos_2", "ref_2", "alt_2", "epcr"]


def run(args: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if check and p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(args)}\n{p.stderr.strip()}")
    return p


def read_contigs(vcf: Path) -> set[str]:
    p = run(["bcftools", "view", "-h", str(vcf)])
    contigs = set()
    for line in p.stdout.splitlines():
        if line.startswith("##contig=<ID="):
            m = re.match(r"##contig=<ID=([^,>]+)", line)
            if m:
                contigs.add(m.group(1))
    if not contigs:
        # Header may omit contigs. Query the index if present.
        p2 = run(["bcftools", "index", "-s", str(vcf)], check=False)
        contigs.update(x.split("\t", 1)[0] for x in p2.stdout.splitlines() if x.strip())
    return contigs


def contig_aliases(chrom: str, contigs: set[str]) -> str | None:
    chrom = chrom.strip()
    options = [chrom]
    if chrom.lower().startswith("chr"):
        options.append(chrom[3:])
    else:
        options.append("chr" + chrom)
    # Support the common mitochondrial aliases without guessing other contigs.
    if chrom in {"M", "MT", "chrM", "chrMT"}:
        options.extend(["M", "MT", "chrM", "chrMT"])
    for c in options:
        if c in contigs:
            return c
    return None


def parse_sample_fields(fields: list[str]) -> list[dict[str, str]]:
    """Parse repeated groups of SAMPLE, GT, DP, GQ, AD after the first four VCF fields."""
    rest = fields[4:]
    samples = []
    for i in range(0, len(rest), 5):
        chunk = rest[i:i+5]
        if len(chunk) < 5:
            continue
        samples.append(dict(zip(["sample", "gt", "dp", "gq", "ad"], chunk)))
    return samples


def genotype_alt_state(gt: str, alt_index: int) -> str:
    if not gt or gt in {".", "./.", ".|."}:
        return "MISSING_GT"
    alleles = re.split(r"[/|]", gt)
    if any(a == "." for a in alleles):
        return "PARTIAL_GT"
    try:
        ints = [int(a) for a in alleles]
    except ValueError:
        return "UNPARSEABLE_GT"
    if alt_index in ints:
        return "ALT_CALLED"
    if all(a == 0 for a in ints):
        return "REF_ONLY"
    return "OTHER_ALT_CALLED"


def ad_for_alt(ad: str, alt_index: int) -> tuple[int | None, int | None]:
    if not ad or ad == ".":
        return None, None
    try:
        values = [int(x) for x in ad.split(",")]
        ref_depth = values[0] if values else None
        alt_depth = values[alt_index] if alt_index < len(values) else None
        return ref_depth, alt_depth
    except ValueError:
        return None, None


def query_locus(vcf: Path, contig: str, pos: int) -> list[list[str]]:
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE\\t%GT\\t%DP\\t%GQ\\t%AD]\\n"
    p = run(["bcftools", "query", "-r", f"{contig}:{pos}-{pos}", "-f", fmt, str(vcf)])
    return [line.split("\t") for line in p.stdout.splitlines() if line.strip()]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, help="Authorized local original VCF/VCF.GZ")
    ap.add_argument("--csv", required=True, help="Track 1 predictions CSV")
    ap.add_argument("--out", required=True, help="Local TSV report path; do not commit to GitHub")
    ap.add_argument("--sample", default=None, help="Optional exact sample ID in the VCF")
    args = ap.parse_args()
    vcf = Path(args.vcf).expanduser().resolve()
    csv_path = Path(args.csv).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    if not vcf.is_file():
        print(f"ERROR: VCF not found: {vcf}", file=sys.stderr); return 2
    if not csv_path.is_file():
        print(f"ERROR: CSV not found: {csv_path}", file=sys.stderr); return 2
    if not shutil.which("bcftools"):
        print("ERROR: bcftools not found on PATH", file=sys.stderr); return 2
    if not (Path(str(vcf)+".tbi").exists() or Path(str(vcf)+".csi").exists() or Path(str(vcf)+".idx").exists()):
        print("ERROR: No adjacent VCF index (.tbi/.csi/.idx). Run `bcftools index -c <VCF>` locally, then rerun.", file=sys.stderr)
        return 2
    try:
        with csv_path.open(encoding="utf-8-sig", newline="") as f:
            r = csv.DictReader(f); rows = list(r); headers = r.fieldnames or []
        missing = [x for x in REQUIRED if x not in headers]
        if missing:
            print("ERROR: CSV missing required fields: " + ", ".join(missing), file=sys.stderr); return 2
        contigs = read_contigs(vcf)
        output = []
        for rowno, row in enumerate(rows, start=2):
            for allele_n in (1, 2):
                chrom = (row.get(f"chrom_{allele_n}") or "").strip()
                pos_raw = (row.get(f"pos_{allele_n}") or "").strip()
                ref = (row.get(f"ref_{allele_n}") or "").strip().upper()
                alt = (row.get(f"alt_{allele_n}") or "").strip().upper()
                if not any([chrom, pos_raw, ref, alt]):
                    continue
                item: dict[str, Any] = {
                    "csv_row": rowno, "proband_id": row.get("proband_id", ""), "allele_number": allele_n,
                    "submitted_chrom": chrom, "submitted_pos": pos_raw, "submitted_ref": ref, "submitted_alt": alt,
                    "vcf_contig": "", "exact_record_and_allele": "NO", "genotype_state": "NOT_CHECKED",
                    "gt": "", "dp": "", "gq": "", "ad": "", "alt_depth": "", "observed_alt_fraction": "",
                    "warnings": ""
                }
                try:
                    pos = int(pos_raw)
                    if pos <= 0: raise ValueError
                except ValueError:
                    item["warnings"] = "INVALID_POSITION"; output.append(item); continue
                vcf_contig = contig_aliases(chrom, contigs)
                if not vcf_contig:
                    item["warnings"] = "NO_MATCHING_CONTIG_NAME"; output.append(item); continue
                item["vcf_contig"] = vcf_contig
                recs = query_locus(vcf, vcf_contig, pos)
                matches = []
                for fields in recs:
                    if len(fields) < 4: continue
                    try: rec_pos = int(fields[1])
                    except ValueError: continue
                    alts = fields[3].upper().split(",")
                    if rec_pos == pos and fields[2].upper() == ref and alt in alts:
                        alt_index = alts.index(alt) + 1
                        matches.append((fields, alt_index))
                if not matches:
                    item["warnings"] = "EXACT_REF_ALT_NOT_FOUND_AT_POSITION"; output.append(item); continue
                fields, alt_index = matches[0]
                item["exact_record_and_allele"] = "YES"
                sample_records = parse_sample_fields(fields)
                if args.sample:
                    sample_records = [x for x in sample_records if x.get("sample") == args.sample]
                elif len(sample_records) > 1:
                    item["warnings"] = "MULTIPLE_SAMPLES; specify --sample for genotype QC"
                    output.append(item); continue
                if not sample_records:
                    item["warnings"] = "ALLELE_FOUND_BUT_SAMPLE_NOT_RESOLVED"; output.append(item); continue
                s = sample_records[0]
                item["gt"], item["dp"], item["gq"], item["ad"] = s["gt"], s["dp"], s["gq"], s["ad"]
                item["genotype_state"] = genotype_alt_state(s["gt"], alt_index)
                ref_depth, alt_depth = ad_for_alt(s["ad"], alt_index)
                item["alt_depth"] = "" if alt_depth is None else alt_depth
                if ref_depth is not None and alt_depth is not None and (ref_depth + alt_depth) > 0:
                    item["observed_alt_fraction"] = f"{alt_depth/(ref_depth+alt_depth):.5f}"
                if item["genotype_state"] != "ALT_CALLED":
                    item["warnings"] = "ALLELE_IN_RECORD_BUT_NOT_CALLED_IN_SELECTED_GENOTYPE"
                output.append(item)
        out.parent.mkdir(parents=True, exist_ok=True)
        cols = ["csv_row", "proband_id", "allele_number", "submitted_chrom", "submitted_pos", "submitted_ref", "submitted_alt",
                "vcf_contig", "exact_record_and_allele", "genotype_state", "gt", "dp", "gq", "ad", "alt_depth", "observed_alt_fraction", "warnings"]
        with out.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(output)
        n = len(output); exact = sum(x["exact_record_and_allele"] == "YES" for x in output)
        called = sum(x["genotype_state"] == "ALT_CALLED" for x in output)
        print(f"Source gate complete: {exact}/{n} submitted alleles match exact VCF record alleles; {called}/{n} are called in the resolved genotype.")
        print(f"Local evidence audit written to: {out}")
        if exact != n or called != n:
            print("WARNING: inspect the local audit and resolve mismatches before upload. Do not publish this TSV.")
            return 1
        print("PASS: exact tuple and genotype-presence checks passed. This is not a pathogenicity or phasing validation.")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr); return 2

if __name__ == "__main__":
    raise SystemExit(main())
