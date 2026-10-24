#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="/mnt/c/MVA_Hackathon/processed"
VCF="/mnt/c/MVA_Hackathon/data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
PAIRS="$BASE/variantforge_v0_5_3/variantforge_v053_pair_candidates.tsv"
HPO="$BASE/variantforge_v0_4_1/gene_disease_summary.tsv"
INHERITANCE="$BASE/variantforge_v0_5/inheritance_overrides_curated.tsv"
OUT="$BASE/pivot_mva_v1_1"

for f in "$VCF" "$PAIRS" "$HPO" "$INHERITANCE"; do
  if [[ ! -f "$f" ]]; then echo "ERROR: required local input missing: $f" >&2; exit 2; fi
done
command -v bcftools >/dev/null || { echo 'ERROR: bcftools not found on PATH' >&2; exit 2; }

python3 "$HERE/src/pivot_mva.py" \
  --pairs "$PAIRS" \
  --hpo-summary "$HPO" \
  --inheritance "$INHERITANCE" \
  --out "$OUT"

CANDIDATE="$OUT/pivot_mva_submission_primary_only.csv"
python3 "$HERE/src/validate_track1_submission.py" "$CANDIDATE"

GATE="$OUT/pivot_mva_source_gate_LOCAL_ONLY.tsv"
python3 "$HERE/src/pivot_mva_source_gate.py" \
  --vcf "$VCF" \
  --csv "$CANDIDATE" \
  --out "$GATE"

printf '\nPreflight finished. Review before upload:\n  %s\n  %s\n  %s\n' \
  "$OUT/RUN_REPORT.txt" "$CANDIDATE" "$GATE"
printf 'Do not commit the local output directory or source-gate TSV to GitHub.\n'
