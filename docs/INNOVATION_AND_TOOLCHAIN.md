# PIVOT-MVA — Innovation and Toolchain Brief

## Innovation claim

**Proof-carrying prioritization:** each nominated pair is accompanied by a structured record of evidence type, evidence provenance, unresolved uncertainty and the next observation most likely to change the ranking. PIVOT-MVA combines pair-aware inheritance reasoning with phenotype-domain counterfactuals and a local exact-allele source gate. Its goal is not merely to predict a gene, but to expose why a pair is ranked, what remains unknown and how the hypothesis can be falsified.

## Four innovation mechanisms

1. **Missingness semantics:** AF not present in a resource remains unknown; no ClinVar match is not benign; unphased heterozygotes are not `in trans`.
2. **Phenotype-domain robustness:** remove one phenotype domain at a time and quantify candidate rank movement; correlated growth or perinatal terms cannot dominate by simple term count.
3. **Evidence-conflict ledger:** separate variant-level assertion, condition association, inheritance model, predicted molecular effect and phenotype fit. Complex susceptibility cannot silently support a Mendelian pair.
4. **Local Source Concordance Gate:** check exact candidate allele tuples and genotype quality fields against the original VCF without uploading a genome-wide table to external services.

## Tool matrix

| Tool | Contribution | Status in this package | Data boundary |
|---|---|---|---|
| Python stdlib | Reproducible ranking, testing, file outputs | Used | Local |
| VariantForge candidate tables | Existing variant annotation/pair discovery | Prior local output; source files in repo | Local TSVs remain private |
| `bcftools` | Exact VCF allele and genotype audit | Source-gate code supplied; user must execute locally | Local VCF only |
| HPO / local HPO association table | Gene–phenotype evidence | Used in prior local workflow | Local/reference files |
| LIRICAL | Independent phenotype-only comparator | Installed locally; output is not EPCR | HPO terms only if external service is avoided |
| Phen2Gene | Orthogonal phenotype-based comparator | Optional targeted check; not run for this artifact | HPO terms only; check current service terms |
| SpliceAI Lookup / Pangolin | Targeted splice-effect review of non-coding candidates | Optional rescue extension; not run for primary pair | Candidate-level lookups only; no whole VCF upload |
| VariantValidator | HGVS/transcript representation checks | Optional pre-submission verification | Minimal candidate-level request |
| ClinVar / gnomAD | Exact allele assertions / reference population frequencies | Local annotations used previously; data coverage is partial | Do not interpret missing as zero |

## Why not train a new neural network?

One proband is not a training cohort. A classifier trained or tuned on one known case would overfit and would not provide meaningful independent validation. The more defensible innovation is a transparent evidence protocol with counterfactual testing, source concordance and explicit uncertainty.

## Judge-facing patient-centered value

The evidence card deliberately distinguishes (a) what supports the hypothesis, (b) what remains unknown, and (c) which additional observation could resolve uncertainty. It avoids converting the hackathon's rank score into a clinical probability or diagnosis. The family/proband distinction is preserved for HPO data and phase is not claimed without segregation evidence.
