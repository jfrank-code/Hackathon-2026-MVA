# PIVOT-MVA v1.0 — Methods and validation status

## Objective

PIVOT-MVA prioritizes candidate variant pairs for one rare-disease proband using existing local annotations. It is a deterministic, evidence-decomposed ranker; it does not predict a calibrated posterior probability and does not diagnose the proband.

## Input and provenance

The method consumes VariantForge's `variantforge_v053_pair_candidates.tsv`, the HPO gene–disease summary generated from the local HPO release, and a local reviewed inheritance table. It never opens or retransmits the raw VCF. The disease-level override file is written separately to preserve provenance and avoid silently assuming that every gene–disease annotation is monogenic.

The HPO information-content constants are copied from the user's pre-existing `hpo_rank_summary.txt` output: rhabdomyosarcoma 6.2328; nephrocalcinosis 4.5080; short stature 1.6942; failure to thrive 2.1723; skeletal muscle atrophy 2.5125; premature birth 4.1259; small for gestational age 3.5359. These are relative information-content values from that local analysis, not disease probabilities.

## Feature architecture

### 1. Association gate

Reviewed disease-level associations are categorized as (a) reviewed Mendelian disease, (b) unclassified/requires review, or (c) complex susceptibility. The default reviewed rules mark OMIM:257300 (MVA1/BUB1B) as a Mendelian disease association with an AR mechanism, sourced to ClinGen; OMIM:212750 is marked complex susceptibility, sourced to GeneReviews. Unknown disease IDs stay unknown. The rules are inspectable in `pivot_mva_association_overrides.tsv` and can be edited with sources.

### 2. Pair/inheritance compatibility

For a reviewed AR model, two heterozygous calls can be compatible with a biallelic hypothesis, but their phase remains unknown unless phase evidence exists in input. The ranker never upgrades two unphased 0/1 calls to confirmed compound heterozygosity. Dominant or unknown inheritance models are not automatically interpreted as a recessive pair.

### 3. Exact clinical evidence

ClinVar P/LP counts only if the per-allele exact-match field supports the match. Disease-label linkage is reported separately from the global assertion. Exact benign and conflicting assertions are separate penalties. No ClinVar match is not assumed benign.

### 4. Phenotype-domain support

The matched HPO profile from the best gene–disease association is parsed into five domains: oncologic, renal, growth, neuromuscular and perinatal. Within a domain only the strongest matched term contributes, reducing double counting among correlated growth/development terms. Exact matches use factor 1.0; specific-subtype matches 0.65; broad-ancestor matches 0.25. The score supports ranking and is not a likelihood ratio. `COMPLEX_SUSCEPTIBILITY` associations do not contribute Mendelian-pair phenotype support.

### 5. Population and technical evidence

AF is used only when an exact allele frequency is supplied; unavailable AF is missing data, not zero. DP, GQ, read allele fraction, nearby indels and MHC location are review flags. Nearby positions do not prove a complex event; MHC location alone does not prove mapping error.

### 6. Ranking policy

Candidates are ordered lexicographically, rather than by an opaque weighted sum:
1. reviewed Mendelian association before unclassified association before complex susceptibility;
2. reviewed inheritance/genotype compatibility;
3. exact benign/conflicting findings;
4. exact ClinVar P/LP and condition linkage;
5. HPO-domain support and domain coverage;
6. LoF/consequence evidence, exact allele rarity and technical quality;
7. nearby-indel/MHC audit flags and stable string keys.

This is an ordinal policy, not a probability model. The ranking and all component axes are written to TSV for audit.

### 7. Counterfactual sensitivity

For the top 50 baseline candidates, the script re-ranks the complete candidate set after removing each phenotype domain in turn. This reveals sensitivity to missing/less-informative domains. It is not a measure of causal probability; it is a robustness diagnostic.

## Validation

The script was syntax-checked and exercised on synthetic fixtures for: BUB1B-like reviewed AR pair, HLA-like complex susceptibility pair, missing AF, ClinVar evidence, phenotype-domain collapse, CSV output, and top-rank behavior. It has not been run against the user's WSL files from this environment; the user must execute it locally and review `RUN_REPORT.txt`.

## Expected current leading hypothesis and boundaries

Prior local VariantForge output already ranked BUB1B first, with a truncating allele `chr15:40209701 T>G` (p.Leu737Ter; exact ClinVar P/LP linked to MVA1) and a rare missense allele `chr15:40220612 T>G` (p.Asn1002Lys; no exact ClinVar match in the local resource). PIVOT-MVA does not encode these coordinates in its source and must reproduce the first rank from the input evidence and configured disease-level rule. This remains a hypothesis until the second allele's effect and phase are independently resolved.

Public participant repositories have also reported the same pair and one reports a full score under the official evaluator. Those reports are external corroboration and are not hard-coded scoring features; their claims do not validate this script or remove the need for manual review.

## EPCR limitation

No independent set of solved cases was used to calibrate EPCR. The primary-only file's default 0.85 and the diverse draft's decreasing secondary values are documented judgement/rank estimates, not trained probabilities. A proper calibration would require independent labeled cases and a held-out calibration/validation strategy. The official evaluator relies on rank/threshold behaviour; the ranking must therefore be reviewed before upload.

## Data governance

The script writes derived individual variant/pair tables only to local output. Do not commit raw/derived variant-level tables or clinical data to a public Git repository. Include only code and public-safe methods documents, and follow the organizer's data-use, retention and deletion instructions.

## References

- Official challenge evaluator: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py
- ClinGen BUB1B/MVA1: https://search.clinicalgenome.org/kb/gene-validity/CGGV%3Aassertion_59147f27-d5a3-4760-ba8d-0429bae3c906-2019-11-22T14%3A53%3A26.352Z
- HPO: https://hpo.jax.org/
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- gnomAD: https://gnomad.broadinstitute.org/
- LIRICAL: https://github.com/TheJacksonLaboratory/LIRICAL
- Exomiser: https://github.com/exomiser/Exomiser
- Public competitor methods: https://github.com/bigbag/mva-hackathon-2026 ; https://github.com/bralewild/mva-hackathon-2026
