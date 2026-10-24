#!/usr/bin/env bash
# Optional, separate validation stage for a small set of already-prioritized candidates.
# It records external predictor outputs; it intentionally does NOT modify PIVOT-MVA scores or EPCR.
set -u
RAW="${RAW_VCF:-/mnt/c/MVA_Hackathon/data/WGS_EX2312012_HGWCNDSX7.vcf.gz}"
OUT="${OUT_DIR:-/mnt/c/MVA_Hackathon/processed/pivot_mva_v1_1/extra_validation}"
mkdir -p "$OUT"

for cmd in bcftools curl python3; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: missing dependency: $cmd" >&2
    exit 2
  fi
done
if [ ! -r "$RAW" ]; then
  echo "ERROR: VCF not readable: $RAW" >&2
  echo "Set RAW_VCF=/path/to/original.vcf.gz if your path differs." >&2
  exit 2
fi

# First check the specific proposed intronic tuple; if it is not present in this exact representation,
# record that fact and do not query splice predictors for it. Fall back to the submitted missense allele.
INTRONIC="$(bcftools query -r chr15:40216470-40216470 -f '%CHROM-%POS-%REF-%ALT\n' "$RAW" 2>/dev/null | grep -Fx 'chr15-40216470-A-G' || true)"
if [ -n "$INTRONIC" ]; then
  VARIANT='chr15-40216470-A-G'
  printf 'Intronic candidate exact tuple found in VCF query: %s\n' "$VARIANT" | tee "$OUT/variant_selection.txt"
else
  VARIANT='chr15-40220612-T-G'
  printf 'Intronic candidate chr15-40216470-A-G not confirmed by exact query; checking submitted missense: %s\n' "$VARIANT" | tee "$OUT/variant_selection.txt"
fi

for TOOL in spliceai pangolin; do
  FILE="$OUT/${TOOL}.json"
  URL="https://${TOOL}-38-xwkwwwxdwq-uc.a.run.app/${TOOL}/"
  echo "Querying ${TOOL} for ${VARIANT} ..."
  if curl -fsS --max-time 20 -G "$URL" \
      --data-urlencode 'hg=38' \
      --data-urlencode "variant=$VARIANT" \
      -o "$FILE"; then
    python3 -m json.tool "$FILE" > "$OUT/${TOOL}.pretty.json" 2>/dev/null || true
    echo "Saved: $FILE"
  else
    echo "WARNING: ${TOOL} query failed; no score inferred." | tee -a "$OUT/external_validation_warnings.txt"
  fi
  sleep 4
done

PHENO="$OUT/phen2gene.json"
HPO='HP:0002859;HP:0000121;HP:0004322;HP:0001508;HP:0003202;HP:0001622;HP:0001518'
echo 'Querying Phen2Gene using the seven proband HPO terms ...'
if curl -fsS --max-time 20 -G 'https://phen2gene.wglab.org/api' \
    --data-urlencode "HPO_list=$HPO" \
    --data-urlencode 'weight_model=sk' \
    -o "$PHENO"; then
  python3 - "$PHENO" <<'PY'
import json, sys
with open(sys.argv[1], encoding='utf-8') as f:
    data = json.load(f)
print('Phen2Gene errors:', data.get('errors', []))
for item in data.get('results', [])[:15]:
    print(item.get('Rank'), item.get('Gene'), item.get('Score'))
for item in data.get('results', []):
    if item.get('Gene') == 'BUB1B':
        print('BUB1B:', item)
        break
PY
else
  echo 'WARNING: Phen2Gene query failed; no ranking inferred.' | tee -a "$OUT/external_validation_warnings.txt"
fi

cat > "$OUT/README.txt" <<EOF
Optional external validation run
Variant queried by SpliceAI/Pangolin: $VARIANT
Original VCF: $RAW
Generated: $(date -Is 2>/dev/null || date)
This stage is an audit only. It does not modify the PIVOT-MVA ranking, pair score, or EPCR.
Keep these local files private; do not commit VCFs, genotype tables, or sample-level audit outputs.
EOF

echo "External validation finished. Outputs: $OUT"
