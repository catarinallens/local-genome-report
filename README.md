# Genome Lab

A private, local tool for exploring your own DNA against public research. It
compares your consumer DNA data (raw genotype files and low-pass whole-genome
sequencing) with the GWAS Catalog, ClinVar, the PGS Catalog and new papers.
It turns the results into a plain-language report you can rerun whenever new
research comes out.

**Your DNA never leaves your computer.** Only public reference databases are
downloaded; all analysis happens locally.

> **Not medical advice.** This is a research and curiosity tool. It can't
> diagnose anything or rule anything in or out. Imputed genotypes are
> statistical estimates, and rare "pathogenic" findings need clinical
> confirmation before they mean anything.

![At a glance](docs/at-a-glance.png)
<sub>Example report built from synthetic test data.</sub>

## What it does

- **Reads your data:**
  - raw genotype files from MyHeritage, 23andMe and AncestryDNA
  - low-pass whole-genome sequencing (CRAM/BAM), which it imputes to about 27
    million variants with [GLIMPSE2](https://github.com/odelaneau/GLIMPSE)
    and the 1000 Genomes reference panel
- **Builds a report in ten parts**: women's & reproductive health, nutrition,
  body/heart/lungs, mind & behaviour, nervous system, hormones & immune
  system, eyes/skin/blood, medicines, ancestry & family, and rare variants &
  carrier status.
- **Labels every card with how established the science is**: strong, moderate
  or early.
- **Single-gene traits** with solid evidence, for example lactase persistence,
  coeliac HLA types, HFE, Lp(a), APOE, alpha-1 antitrypsin, CFH/ARMS2, and
  bitter taste.
- **About 30 published polygenic scores** from the PGS Catalog, shown as an
  estimated lifetime chance next to the typical one, with percentiles among
  1000 Genomes reference people.
- **Medicine response (pharmacogenomics)** for CYP2C19, CYP2C9, VKORC1,
  SLCO1B1, DPYD, TPMT, NUDT15, CYP3A5, UGT1A1, CYP2B6, partial CYP2D6 and more,
  with what CPIC guidelines say, including for antidepressants.
- **ClinVar findings**, split into "could matter for you" and "what you could
  pass on".
- **Ancestry check** (PCA against 1000 Genomes) and **maternal haplogroup**
  (Haplogrep 3), when your file contains mitochondrial DNA.
- **Watch list**: follow conditions you care about, with new studies and
  papers flagged in a "What's new" section each time you rerun it.

![Nutrition](docs/nutrition.png)
![Medicines](docs/meds.png)

## Requirements

macOS or Linux, Python 3.7+, and a few bioinformatics tools from
[conda](https://docs.conda.io/) / bioconda:

    conda create -y -n genome -c conda-forge -c bioconda python samtools bcftools glimpse-bio
    conda install -y -n genome -c conda-forge -c bioconda plink2 openjdk=17 minimap2
    conda activate genome

Disk space: about 13 GB once everything is set up, and up to about 15 GB while
imputing. Data lives in `~/genome-lab`; set `GENOME_LAB_HOME` to change that.

## Quick start

With a raw genotype file only (takes minutes):

    python3 genome_lab.py sex female            # or male: decides X-chromosome handling and which studies are preferred
    python3 genome_lab.py ancestry european     # or "any"
    python3 genome_lab.py chip YOUR_RAW_DATA.zip
    python3 genome_lab.py run                   # downloads the databases and opens the report

With low-pass sequencing (CRAM), add:

    python3 genome_lab.py setup                                   # reference genome, ~900 MB
    caffeinate -i python3 genome_lab.py impute YOUR_FILE.cram     # 6–12 hours; resumable
    caffeinate -i python3 genome_lab.py run                       # scores, ancestry, medicines, report

(`caffeinate -i` keeps a Mac awake; leave it out on Linux.) Every long step
can be stopped with Ctrl+C and resumed by running the same command again.

## Commands

| Command | What it does |
|---|---|
| `run` | Refresh databases (weekly), compute what's missing, build a new report |
| `report` | Rebuild the report only |
| `lookup rs7903146` | Your genotype at a variant, plus every study and ClinVar entry about it |
| `watch "Migraine" migraine` | Follow a condition (studies + papers) |
| `unwatch "Migraine"` | Stop following it |
| `score "Breast cancer" PGS000004` | Choose a specific published score (`auto` or `none` also work) |
| `extras --force` | Recompute scores, medicines and traits |
| `maternal FILE` | Read the maternal haplogroup from a file that contains mtDNA |
| `status` | What's done so far |

## How it works, briefly

- **Positions** are matched on GRCh37. Results from the GWAS Catalog (GRCh38)
  are lifted over with the Ensembl chain file.
- **Catalog variant counts** keep one signal per region (±500 kb clumping).
  Effect sizes are adjusted by how precisely each study measured them
  (empirical Bayes shrinkage). These counts are shown only as "more or fewer
  than typical"; published scores are used for anything resembling risk.
- **Polygenic scores** use the PGS Catalog's harmonised GRCh37 files and
  [plink2](https://www.cog-genomics.org/plink/2.0/). Percentiles come from 1000
  Genomes reference people, as in
  [pgsc_calc](https://github.com/PGScatalog/pgsc_calc). Lifetime chances
  combine population prevalence for women of European ancestry with each
  score's effect size, so they're approximate.
- **Medicine phenotypes** follow CPIC rules from single variants. Gene
  deletions, duplications and small insertions/deletions (for example CYP2D6
  \*3, \*5 and \*6) aren't detected.

## Validation

- **Imputation**: on held-out 1000 Genomes people, genotype concordance was
  about 98–99% for common variants (chromosome 22 and X). Polygenic scores
  from imputed data were within about 0.1 SD of scores from true genotypes.
- **Single-gene and medicine calls**: run on all 2,504 1000 Genomes people,
  they reproduce published population frequencies (for example lactase
  persistence, HERC2 eye colour, ADH1B, and CPIC phenotype frequencies for
  CYP2C19, CYP2C9, SLCO1B1, TPMT, DPYD).
- **Raw-file formats**: 23andMe, AncestryDNA and MyHeritage formats are
  tested, including no-calls, X and mitochondrial lines.

## Limitations

- Common-variant results shift risk by small amounts and capture only part of
  each condition's genetics. They're far less informative than symptoms and
  family history.
- Very rare variants are poorly imputed from low-pass data; a "not found" is
  not a negative test.
- Most reference data and scores come from European-ancestry studies, so
  results are less accurate for other ancestries.
- CNVs, repeat expansions, most indels and HLA typing (estimated here from tag
  variants) need clinical tests.

## Data sources

[GWAS Catalog](https://www.ebi.ac.uk/gwas/) ·
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) ·
[PGS Catalog](https://www.pgscatalog.org/) ·
[1000 Genomes](https://www.internationalgenome.org/) ·
[Europe PMC](https://europepmc.org/) ·
[CPIC](https://cpicpgx.org/) ·
[Haplogrep 3](https://github.com/genepi/haplogrep3) /
[PhyloTree](https://www.phylotree.org/)

All of them are downloaded directly from their providers on your machine;
none are redistributed here.

## License

MIT. See [LICENSE](LICENSE).
