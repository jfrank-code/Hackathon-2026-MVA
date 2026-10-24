# PIVOT-MVA — Track 1 Methods Report

**Challenge:** Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1 (Variant Prediction)  
**Method:** PIVOT-MVA, with a source-concordance gate and counterfactual evidence audit  
**Output:** one primary compound-heterozygous hypothesis  
**Genome build:** GRCh38, as recorded in the existing local annotation workflow

## Executive summary

PIVOT-MVA is a lightweight, evidence-decomposed prioritization approach for a single-proband rare-disease case. Rather than treating every annotation as an interchangeable score, it separates technical call quality, exact-allele clinical evidence, population frequency, gene–disease mechanism, genotype/inheritance compatibility and phenotype concordance. It also records what remains unknown and tests ranking sensitivity when a phenotype domain is removed.

The current primary submission is the pair `chr15:40209701 T>G` and `chr15:40220612 T>G` in **BUB1B**, selected from prior local VariantForge outputs and their evidence review. The first variant was previously annotated as `NM_001211.6:c.2210T>G, p.Leu737Ter` with an exact-allele ClinVar P/LP assertion associated with MVA1. The second was annotated as `NM_001211.6:c.3006T>G, p.Asn1002Lys`; it is rare in the local population-frequency annotation but has no exact ClinVar pathogenic assertion in that resource. These facts make the pair a prioritized hypothesis, not proof that both variants are in trans or independently functionally validated.

**Validation status from the local run:** PIVOT-MVA v1.0 was executed against the locally generated pair-candidate TSV, HPO summary and curated inheritance table. It scored 23,485 pair hypotheses and read 3,845 gene–disease/HPO associations; BUB1B ranked first. The Source Concordance Gate then checked the two submitted alleles against the original local VCF: **2/2 exact CHROM/POS/REF/ALT tuples matched VCF records and 2/2 were called in the resolved sample genotype**. This verifies allele-record concordance and genotype presence only; it does not verify pathogenicity, independent clinical causality, transcript consequence, calibration of EPCR, or phase. PIVOT-MVA uses the prior annotation/pair tables for ranking; it does not re-annotate every base in the raw VCF itself.

**Additional targeted external checks (run 2026-10-24):** the exact tuple `chr15:40216470 A>G` proposed as a deep-intronic candidate was not confirmed by the exact-position/allele lookup against this VCF query and was not scored. The submitted missense allele `chr15:40220612 T>G` was queried in the GRCh38 SpliceAI Lookup API. For the primary BUB1B transcript `ENST00000287598.11`, the SpliceAI delta scores were `DS_AG=0.016`, `DS_AL=0.001`, `DS_DG=0.003`, `DS_DL=0.000` (maximum 0.016); Pangolin returned `DS_SG=0.010` and `DS_SL=-0.003`. These are low predicted splicing effects for this allele, not evidence of protein-level benignity or pathogenicity. Phen2Gene was run on the seven HPO terms as a phenotype-only comparison; BUB1B ranked 26 (`Score=0.435167`), while NBN, PTH1R and KRAS ranked 1–3. Thus, this independent phenotype-only tool did **not** strongly prioritize BUB1B; this discordance is reported rather than hidden. These API results were not used to modify PIVOT-MVA's ranking or the submitted EPCR.

## 1. Data and phenotype representation

The starting genomic input is the challenge-provided WGS VCF (approximately 5.0 million records in the previous local inventory). Candidate annotation and pair generation were performed in the user's local workflow, with the MANE-annotated VCF and local HPO/ClinVar/population references. The submitted repository does not include the VCF, clinical source document, genome-wide genotype tables, or private intermediate files.

The seven proband phenotype terms used in the existing local workflow are:

- Rhabdomyosarcoma (`HP:0002859`)
- Nephrocalcinosis (`HP:0000121`)
- Short stature (`HP:0004322`)
- Failure to thrive (`HP:0001508`)
- Skeletal muscle atrophy (`HP:0003202`)
- Premature birth (`HP:0001622`)
- Small for gestational age (`HP:0001518`)

The parents' recurrent spontaneous abortions (`HP:0200067`) are retained as family-history context, not incorrectly assigned as a phenotype of the proband.

## 2. Evidence architecture

### 2.1 Variant-level quality

Where present in the input annotations, evaluate `FILTER`, `GT`, `DP`, `GQ`, `AD` and allele fraction. A high genotype quality is evidence about the caller's confidence, not proof that a variant is biologically causal. An allele that is present in the VCF is not accepted as a match unless the submitted chromosome, position, reference and alternate allele agree exactly.

### 2.2 Exact-allele clinical evidence

ClinVar assertions are matched to the exact allele rather than inferred from the gene name. Global classification, condition link, review status, conflicting assertions and missing matches are kept distinct. Absence of an exact ClinVar record is not treated as benign evidence.

### 2.3 Population-frequency missingness

Allele frequency is used only when the local annotation contains a relevant match. A missing frequency is “unknown/not observed in this resource,” not AF=0. The existing gnomAD resource was only partial (chr15 exomes noted locally); it is not described as genome-wide coverage.

### 2.4 Condition-specific inheritance and pair state

The pair is compatible with a recessive BUB1B/MVA1 hypothesis because both variants were previously called heterozygous and the local curated gene–condition model is autosomal recessive. However, the two calls do **not** by themselves establish a compound-heterozygous state: phase and parental segregation were not available in the local callset. PIVOT-MVA leaves phase `UNKNOWN` rather than upgrading the pair to confirmed `in trans`.

### 2.5 Phenotype domains instead of naive HPO counting

The PIVOT-MVA prototype collapses related HPO matches into five domains: oncologic, renal, growth, neuromuscular and perinatal. Within a domain, the most informative match contributes rather than allowing correlated terms to inflate the score. A broad ontology relation, a specific-subtype relation and an exact match remain distinguishable. Gene–disease links flagged as complex susceptibility are not allowed to masquerade as Mendelian support for a recessive pair.

### 2.6 Counterfactual robustness

The prototype re-ranks the candidate set while removing one phenotype domain at a time. The resulting rank movement is a robustness diagnostic: it reveals whether a candidate depends on a single phenotype cluster. It is not a causal probability or an independent validation cohort.

## 3. Innovation: evidence-traceable and falsifiable prioritization

The main innovation is a practical combination of four safeguards, rather than a new black-box model trained on one case:

1. **Evidence ledger:** evidence is partitioned by source and type; “unknown,” “negative,” “conflicting” and “not covered” are distinct states. This makes the reasons for a rank auditable.
2. **Counterfactual phenotype audit:** the system tests whether the ranking survives removal of each phenotype domain, making sensitivity to clinical input explicit.
3. **Pair-aware uncertainty:** genotype compatibility, disease mechanism and phase confidence are represented separately. Two heterozygous variants are not called a confirmed biallelic diagnosis without phase/segregation evidence.
4. **Source Concordance Gate:** `pivot_mva_source_gate.py` was run locally with `bcftools`. It returned PASS for exact allele-tuple matching and called-genotype presence for both submitted alleles (2/2 each). The detailed audit contains sample-level fields and must stay local; only the aggregate check result is reported here. This gate does not establish pathogenicity or phase.

This design aims to make a small ranking more defensible and reusable across cases, especially where annotation gaps, repetitive regions, complex-susceptibility associations or correlated symptoms can produce false confidence.

## 4. Targeted external validation beyond the core ranking

### 4.1 SpliceAI Lookup and Pangolin — queried

The deep-intronic tuple `chr15:40216470 A>G` mentioned as a candidate in a public participant report was not confirmed by the local exact-allele check at that coordinate and therefore was not submitted to splice prediction. This result only says that the exact tuple was not found by the specific query; it does not prove absence under all possible reference-build or variant representations.

The submitted BUB1B missense allele `chr15:40220612 T>G` was queried with the GRCh38 APIs. For primary transcript `ENST00000287598.11`, SpliceAI reported `DS_AG=0.016`, `DS_AL=0.001`, `DS_DG=0.003`, and `DS_DL=0.000` (maximum delta score 0.016). Pangolin reported `DS_SG=0.010` and `DS_SL=-0.003`. These low values do not predict a substantial splice effect for this missense allele in the queried transcript. They do **not** evaluate or rule out a protein-level effect of p.Asn1002Lys, and they cannot classify the variant as benign. SpliceAI Lookup provides these model scores for user-specified variants and recommends individual queries rather than high-volume API use ([service documentation](https://spliceailookup.broadinstitute.org/); [API documentation](https://github.com/broadinstitute/SpliceAI-lookup)).

### 4.2 Phen2Gene — phenotype-only contrast, not concordant

Phen2Gene was queried using the same seven proband HPO terms and `weight_model=sk`; the API returned no errors. It ranked `NBN` first (score 1.000000), `PTH1R` second (0.975613), `KRAS` third (0.975235), `HRAS` fourth (0.918504), and `BUB1B` 26th (0.435167). This is a genuine disagreement with PIVOT-MVA's pair-level ranking: Phen2Gene uses HPO terms for gene prioritization and does not incorporate the submitted pair's exact genotype, ClinVar allele evidence, population frequency, or phase. The result therefore does not invalidate a genotype-driven BUB1B hypothesis, but it means the submission must not claim independent phenotype-only convergence on BUB1B. Phen2Gene's score is a gene-ranking score, not EPCR or a clinical probability ([Phen2Gene API documentation](https://phen2gene.wglab.org/docs)).

### 4.3 Interpretation and unrun tools

These external outputs were recorded for transparency but were **not** added to the numerical PIVOT-MVA score or used to alter `epcr=0.85`. The following checks remain unrun for this submission: AlphaMissense/PrimateAI-3D/AlphaGenome for protein-level impact; VariantValidator for independent HGVS validation; any read-backed or parental phase/segregation analysis. No hosted service received the complete VCF or genome-wide genotype table.

## 5. Candidate evidence card

### Supporting observations from prior local outputs

- Two BUB1B variants were previously present as heterozygous calls in the annotation workflow.
- The truncating allele `chr15:40209701 T>G` was reported with the exact-transcript annotation `p.Leu737Ter` and an exact-allele ClinVar P/LP assertion associated with MVA1.
- The second allele `chr15:40220612 T>G` was reported as `p.Asn1002Lys`, with high call-quality values in the prior annotation output and a very low frequency in the partial local gnomAD annotation.
- The local disease association review classified the BUB1B/MVA1 model as autosomal recessive.

### Unresolved evidence

- Phase is unknown; two heterozygous calls do not establish `in trans`.
- Parental genotypes and segregation were not available in the records used for this report.
- No functional assay is reported for the missense allele in this workflow.
- The partial local gnomAD file must not be treated as a genome-wide search.
- The Source Concordance Gate passed the limited exact-allele and genotype-presence checks for both submitted variants. Biological causality and phase remain unresolved.

### Most informative next checks

1. If authorized reads or family genotypes become available, evaluate allele balance and phase/segregation; do not infer `in trans` from unphased `0/1` calls.
2. Recheck the transcript and condition-specific ClinVar evidence for both variants; the exact-record gate confirms coordinates/genotype, not the truth of clinical assertions.
3. For newly prioritized intronic/splice-region candidates, consult SpliceAI/Pangolin and verify transcript/coordinate consistency before including them.

### Falsification ledger: what would change or reject this hypothesis?

- **Source mismatch:** if either submitted `CHROM/POS/REF/ALT` tuple is not an exact allele in the source VCF, stop and correct the representation before submission.
- **Phase contradiction:** if read-backed or parental evidence shows both BUB1B variants are in cis and no other pathogenic allele is found, the proposed recessive pair is weakened substantially.
- **Second-allele contradiction:** if a reliable, condition-relevant functional/clinical review supports a benign interpretation for `p.Asn1002Lys`, do not keep treating rarity as pathogenic evidence.
- **Phenotype fragility:** if the pair's rank collapses when one phenotype domain is removed, describe it as more sensitive to phenotype selection instead of presenting its rank as robust.
- **Reference-coverage gap:** absence from an incomplete local gnomAD file is not evidence that an allele is absent from population databases.

This explicit falsification ledger is a core part of the patient-centered design: it says not only why a pair is prioritized but also what evidence could overturn it.

## 6. Submitted confidence field

The CSV uses `epcr = 0.85` for a single primary row. This is a **judgement-based ordering value, not a calibrated clinical probability**: no independent solved-case cohort or calibration procedure was used. The official scoring method ranks by EPCR and evaluates the candidate variants; the single-row submission avoids diluting precision with weak alternatives. The value must not be communicated to clinicians or families as an 85% diagnostic probability.

## 7. Reproducibility and data governance

The package uses Python's standard library plus local `bcftools` for the source gate; no GPU or newly downloaded large reference database is needed. Existing local annotation/ranking TSV files can be used to reproduce the current PIVOT reranking, while the original VCF remains local. Do not commit VCF/BCF/BAM/CRAM/FASTQ files, indexes, clinical documents, per-variant genotype tables, full candidate evidence tables, local audit outputs or caches.

The official organizer states that Track 1 submissions are eligible for the Innovation Prize. This report therefore describes concrete method contributions, limitations and reproducibility, rather than claiming that the method is a new clinical test. See the official discussion on prizes and recognition: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/28

## 8. Limitations

- Single proband; no supervised model training or independent calibration dataset.
- PIVOT-MVA v1.0 was run on locally generated candidate/HPO/inheritance TSVs, and the selected pair passed the limited source-concordance check against the original VCF. The ranking stage was not a fresh annotation pass across all 5 million raw VCF records.
- The source gate confirms the submitted allele tuples and called genotype, not pathogenicity, phase, or the accuracy of every upstream annotation.
- Phenotype ontology evidence may be incomplete and correlated; domain collapse reduces but does not eliminate this problem.
- Computational predictions and ClinVar records may have limitations or review conflicts; exact variant-level evidence must be interpreted in condition context.
- Independent phenotype-only Phen2Gene ranking placed BUB1B at rank 26 (score 0.435167), so it did not independently support BUB1B as a top phenotype match. This discordance is an uncertainty, not confirmation of or a definitive exclusion of the genotype-based hypothesis.
- This is a hackathon prioritization result, not a clinical diagnosis, independent clinical validation, or treatment recommendation.

## References

1. [Official MVA Hackathon Track 1 evaluator](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)
2. [Official challenge rules, including privacy and data-handling requirements](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)
3. [Organizer statement: Track 1 and Innovation Prize eligibility](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/28)
4. [ClinGen gene–disease validity framework](https://clinicalgenome.org/curation-activities/gene-disease-validity/)
5. [Human Phenotype Ontology](https://hpo.jax.org/)
6. [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
7. [gnomAD](https://gnomad.broadinstitute.org/)
8. [LIRICAL](https://github.com/TheJacksonLaboratory/LIRICAL)
9. [SpliceAI Lookup](https://spliceailookup.broadinstitute.org/)
10. [Phen2Gene](https://phen2gene.wglab.org/)
11. [VariantValidator](https://variantvalidator.org/)
