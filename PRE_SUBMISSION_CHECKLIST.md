# PIVOT-MVA — final submission checklist

## Files for the Track 1 form

- `submission/PIVOT-MVA_track1_predictions.csv`
- `submission/PIVOT-MVA_track1_report.md`
- The real GitHub URL after you create the repository and push the public-safe code.

Suggested Team / Display Name: `PIVOT-MVA — Evidence Ledger`.

## Before uploading

1. Re-run the structural validator:

   ```bash
   python3 src/validate_track1_submission.py submission/PIVOT-MVA_track1_predictions.csv
   ```

2. Re-run the current PIVOT model on the local VariantForge tables and inspect its `RUN_REPORT.txt`.
3. Check the submitted exact allele tuples against the original VCF locally:

   ```bash
   python3 src/pivot_mva_source_gate.py \
     --vcf /path/to/WGS_EX2312012_HGWCNDSX7.vcf.gz \
     --csv submission/PIVOT-MVA_track1_predictions.csv \
     --out /path/to/local-only/pivot_mva_source_gate_LOCAL_ONLY.tsv
   ```

   The VCF needs an index and `bcftools`. If the source-gate TSV reports anything other than exact match plus ALT called, inspect the discrepancy before upload. Never commit that local audit file.

4. The submitted `epcr=0.85` is an uncalibrated ranking estimate. Do not call it an 85% clinical probability.
5. Do not add extra hypotheses to fill ten slots. Additional variants can reduce F-max precision.
6. Do not upload raw VCFs, derived genotype tables, clinical documents or index files to GitHub.
7. After the submission deadline, follow the organizer's requirement to publish the repository. Follow the organizer's data deletion and confirmation instructions within 30 days of close.

## Status that must remain transparent

The previous local VariantForge output had BUB1B as the leading pair. The PIVOT-MVA code is tested on synthetic fixtures, but the current source-gate check and full rerun against the original VCF have not been executed in this artifact-generation environment. Claim neither as passed until you run them locally.
