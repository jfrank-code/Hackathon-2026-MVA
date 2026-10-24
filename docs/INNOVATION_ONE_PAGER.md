# PIVOT-MVA — an audit-first variant-prioritization method

**Problem:** Rare-disease pipelines can produce a plausible ranked gene while hiding why evidence was missing, which annotations conflict, or how much the ranking depends on a correlated group of symptoms.

**Idea:** make a variant pair “carry its own evidence.” PIVOT-MVA separates variant call quality, exact allele-level clinical assertions, population frequency, disease mechanism, inheritance and phenotype fit; a phenotype-domain counterfactual shows how much each domain drives rank. A local source-concordance gate then checks each submitted allele against the original VCF before upload.

**Distinctive contribution:** not another opaque weighted score, and not a neural network trained on one patient. It is a transparent pair-level audit process that makes uncertainty and falsifiers explicit:

- Unknown frequency remains unknown; no exact ClinVar record is not treated as benign.
- Two heterozygous calls do not prove `in trans` phase.
- Complex-susceptibility gene–disease links do not automatically count as Mendelian support.
- The family reproductive-history term remains separate from the proband's own phenotype set.
- Removing one phenotype domain at a time tests robustness rather than rewarding symptom count.
- Source allele tuples and genotype presence are checked locally; sample-level audit outputs never enter the public repository.

**Patient-centered output:** the evidence card communicates what supports a candidate, what remains unresolved and which next observation is most informative. It deliberately avoids turning EPCR into a clinical probability or diagnosis.

**Scalability:** Python standard library, existing local annotations and `bcftools`; targeted optional splice/HPO checks are reserved for candidates that warrant them. No GPU, new multi-gigabyte database or cohort-specific supervised training is required.

**Current evidence ceiling:** BUB1B was the leading pair in prior local VariantForge outputs. The PIVOT-MVA reranking and source-concordance wrapper have synthetic/unit tests, but the packaged run still requires local execution against the authorized source VCF before claiming an end-to-end pass.
