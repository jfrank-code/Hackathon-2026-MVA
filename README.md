# PIVOT-MVA: evidence-traceable variant-pair prioritization

PIVOT-MVA is a lightweight, auditable prioritization workflow for Track 1 of the Rare Disease, Real Kid: MVA Hackathon 2026.

## What is innovative about it?

- Separates exact-allele clinical evidence, population frequency, variant-call quality, phenotype match, inheritance compatibility and phase uncertainty.
- Groups correlated HPO observations into phenotype domains, then tests rank stability by removing domains one at a time.
- Distinguishes missing evidence from negative evidence and flags complex-susceptibility relationships instead of treating them as monogenic support.
- Adds a Source Concordance Gate to compare each submitted `CHROM/POS/REF/ALT` against the original VCF locally using `bcftools`.
- Provides a concise, patient-centered evidence card: what supports the hypothesis, what remains unknown, and which next check could falsify or strengthen it.

This is not a new clinical test, and EPCR is not calibrated. The current submitted candidate came from the existing local VariantForge evidence review. The source gate must be run locally before claiming the exact input genotype check passed.

## Package layout

- `src/pivot_mva.py` — phenotype-domain and pair-aware reranking layer; consumes the existing pair/HPO/inheritance TSVs.
- `src/pivot_mva_source_gate.py` — local validation of submitted alleles against the original VCF with `bcftools`.
- `src/variantforge_v05_pairtrace.py`, `src/variantforge_v053_pairtrace.py`, `src/variantforge_v054_association_aware_rerank.py` — prior pair-generation and association-review layers. They consume locally generated annotation tables; they do not replace upstream annotation of the raw VCF.
- `docs/METHODS_ENGINE.md` — detailed PIVOT-MVA engine methods.
- `src/validate_track1_submission.py` — structural CSV validation.
- `submission/PIVOT-MVA_track1_predictions.csv` — one primary pair for upload.
- `submission/PIVOT-MVA_track1_report.md` — methods report for upload.
- `docs/INNOVATION_AND_TOOLCHAIN.md` — innovation rationale and optional second-pass tools.

## Reproduce the existing local reranking in WSL

Use the data-derived TSVs that exist in your authorized local environment; do not move them into this public repository.

```bash
python3 src/pivot_mva.py \
  --pairs /mnt/c/MVA_Hackathon/processed/variantforge_v0_5_3/variantforge_v053_pair_candidates.tsv \
  --hpo-summary /mnt/c/MVA_Hackathon/processed/variantforge_v0_4_1/gene_disease_summary.tsv \
  --inheritance /mnt/c/MVA_Hackathon/processed/variantforge_v0_5/inheritance_overrides_curated.tsv \
  --out /mnt/c/MVA_Hackathon/processed/pivot_mva_v1_1
```

The path under `/mnt/c` is local and must not be committed. Inspect `RUN_REPORT.txt` and the ranking tables before selecting a CSV. This reranking is not a substitute for preprocessing/annotation of the original VCF.

## Run the source-concordance gate before submission

```bash
python3 src/pivot_mva_source_gate.py \
  --vcf /mnt/c/MVA_Hackathon/data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
  --csv submission/PIVOT-MVA_track1_predictions.csv \
  --out /mnt/c/MVA_Hackathon/processed/pivot_mva_source_gate_LOCAL_ONLY.tsv
```

The VCF must have a `.tbi` or `.csi` index and `bcftools` must be installed. If there is no index, create one locally with `bcftools index -c /path/to/original.vcf.gz`. If a multi-sample VCF is used, add `--sample SAMPLE_ID`. **The output audit contains individual genotype-level fields and is local-only; never commit it.** A PASS verifies source allele identity/genotype presence, not pathogenicity or phase.

## Validate the submission file

```bash
python3 src/validate_track1_submission.py submission/PIVOT-MVA_track1_predictions.csv
```

This validator checks structure only. It does not prove that the candidate is causal or that the source gate passed.

## Optional, targeted tools (not required to run the primary submission)

- [SpliceAI Lookup](https://spliceailookup.broadinstitute.org/) / [Pangolin](https://github.com/tacazares/Pangolin): follow up plausible splice/non-coding variants individually.
- [Phen2Gene](https://phen2gene.wglab.org/) alongside local LIRICAL: compare phenotype-based prioritization from HPO terms.
- [VariantValidator](https://variantvalidator.org/): validate transcript/HGVS representations.
- Local ClinVar and gnomAD: keep exact-match provenance and reference coverage/version in the audit.

Do not upload the full VCF or genotype-level tables to hosted prediction services. These optional lookups have not been used as evidence for the submitted pair unless you independently execute them and document the results.

## Privacy and repository hygiene

The challenge data are gated. Do not commit the original VCF/FASTQ/BAM/CRAM, indices, clinical source documents, genotype-level tables, full per-variant candidate tables or local source-gate output. The organizer has specified deletion of gated/raw and reconstructive data within 30 days of the close. Follow the current organizer instructions.

## Validation status

- PIVOT-MVA prototype: syntax checked and tested on synthetic fixtures.
- Existing local VariantForge result: BUB1B was the top pair in prior local outputs.
- Local PIVOT-MVA v1.0 rerun: **completed** on the existing pair/HPO/inheritance TSVs; 23,485 pair hypotheses scored, 3,845 associations read, BUB1B ranked first. This is a reranking of precomputed candidates, not a fresh annotation pass through every raw VCF record.
- Source-concordance gate: **PASS for the two submitted alleles**; 2/2 exact VCF allele tuples matched and 2/2 were called in the sample genotype. This does not validate pathogenicity or phase. The detailed source-gate TSV and other local outputs must remain outside GitHub.

## License

CC BY 4.0 for this original work, subject to the competition's terms; third-party tools and datasets retain their own licenses.

## Publish the reproducible repository

1. Create an empty repository in GitHub named `pivot-mva-rare-disease` (choose private while the event remains open if you prefer; the organizer says it must be made public after closure).
2. From this package directory, run:

```bash
git init
git add README.md LICENSE .gitignore src tests docs submission
git status --short
# Check that no VCF/BCF/BAM/CRAM/FASTQ, clinical document, or local TSV was staged.
git commit -m "Add PIVOT-MVA evidence-traceable variant prioritization"
git branch -M main
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/pivot-mva-rare-disease.git
git push -u origin main
```

Replace `YOUR_GITHUB_USERNAME` with your real GitHub username and create that empty remote first. The URL is a suggested repository name, not a repository I created. Use the resulting real URL in the required form field. Do not commit local source-gate output or private variant tables.

## One-command local preflight

After unzipping this package inside WSL, execute:

```bash
./run_local_preflight.sh
```

The wrapper reruns the PIVOT layer using your already-generated local TSVs, validates the produced primary-only CSV and checks its allele tuples/genotypes against the original VCF. It intentionally writes all detailed outputs under `/mnt/c/MVA_Hackathon/processed/pivot_mva_v1_1` and does not upload or copy data anywhere. A successful preflight is the point at which you can compare the newly generated candidate with the packaged submission CSV. If they differ, inspect `RUN_REPORT.txt` and the source gate before choosing which CSV to upload.

## Optional external validation stage

After `run_local_preflight.sh`, an optional, separate check can query SpliceAI Lookup, Pangolin and Phen2Gene for the shortlisted variant/phenotype set:

```bash
./scripts/run_external_validation.sh
```

This script requires `bcftools`, `curl`, and `python3`, plus network access. It first checks a proposed intronic tuple against the local VCF and falls back to the submitted missense tuple if that exact representation is not found. It writes JSON outputs under `/mnt/c/MVA_Hackathon/processed/pivot_mva_v1_1/extra_validation` by default. These external results are recorded as a validation audit only; they do not modify PIVOT-MVA ranking or EPCR because no calibration study was performed to combine them numerically. Do not commit the output directory, source VCF, or genotype-level audit files.
