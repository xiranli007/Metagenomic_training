# Metagenomic Analysis Pipeline — Leonard et al. (2025)

**Paper:** "Air microbiomes reveal presence of Shiga toxin-producing *Escherichia coli* in airborne cattle pen soil adjacent to large feedlot"  
**Journal:** Science of the Total Environment 1000 (2025) 180375  
**GitHub:** [mmammel8/kmer_id](https://github.com/mmammel8/kmer_id)

---

## Pipeline Overview

This paper uses shotgun metagenomic sequencing of air and cattle pen soil samples to investigate airborne STEC transfer from a cattle feedlot. The bioinformatic pipeline has **six major modules**:

1. **Read QC & Trimming** — FastQC + Trimmomatic
2. **E. coli Virulence Gene & Serogroup Detection** — BLAST + SeqMan NGen
3. **Bacterial Community Profiling (Taxonomy)** — Custom k-mer classifier (kmer_id)
4. **Mitochondrial DNA Profiling** — Custom k-mer classifier + BURST
5. **Source Tracking** — SourceTracker2
6. **Statistical & Comparative Analysis** — LEfSe, diversity metrics, correlations (R/Python)

---

## Step-by-Step Summary

### Step 1: Read Quality Control & Trimming

**What it does:** Assesses raw read quality and trims adapters/low-quality bases from Illumina paired-end (150 bp) reads generated on NextSeq 550 or NextSeq 2000.

**Tools:**
- FastQC v0.12.1 — quality assessment
- Trimmomatic 0.38 — adapter removal (NexteraPE adapter file) and quality trimming

**Input:** Raw paired-end FASTQ files (R1/R2)  
**Output:** Trimmed paired-end FASTQ files

---

### Step 2: E. coli Virulence Gene & Serogroup Screening

**What it does:** Screens trimmed reads against a custom database of E. coli virulence genes (stx1, stx2, eae, ehxA, saa, subAB, aafA, aggR, bfpA, perA, papD, ipaH, sta, stb, eltAB) and serogroup determinants (wzx, wzy, wzm, wzt, O-group modifier genes). Uses a two-pass approach: (1) low-stringency BLAST to pre-filter candidate reads (≥90% identity over ≥50 nt), then (2) high-stringency filtering and subtyping with SeqMan NGen v18 (DNASTAR, commercial).

**Tools:**
- BLAST+ (blastn) — initial read filtering
- SeqMan NGen v18 (DNASTAR) — high-stringency subtyping (commercial; no open-source equivalent used in the paper)

**Input:** Trimmed FASTQ files + custom virulence/serogroup gene database (FASTA)  
**Output:** Per-sample virulence gene presence/absence table; serogroup assignments

> **Note:** SeqMan NGen is proprietary software. For a fully open-source reproduction, you could substitute the second-pass step with a stricter BLAST or short-read aligner (e.g., BWA-MEM or Bowtie2) against the same database, followed by custom filtering scripts. The Snakemake rules below implement the BLAST pre-filter step; the SeqMan NGen step would be run outside Snakemake or replaced.

---

### Step 3: Bacterial Community Profiling (k-mer Taxonomy)

**What it does:** Assigns taxonomy to trimmed reads using a custom C++ k-mer signature program and a database of 30 bp unique bacterial taxon sequences (~5,900 taxa, avg ~40,000 k-mers each). Low-level identifications are confirmed by BLAST of full-length reads. Subspecies-level relative abundances are rolled up to species level. Taxa with RA ≥ 0.02% are considered present.

**Tools:**
- `nk10` (compiled from `newkmer_10nx.cpp` in [kmer_id repo](https://github.com/mmammel8/kmer_id)) — k-mer based read classification
- `readbatch_10.py` — collects results and generates a CSV report of read counts and % abundance per species
- BLAST+ (blastn) — confirmatory analysis of low-abundance hits

**Database files** (placed in `bact10/` subdirectory):
- `bData10.txt`, `btree10.txt`, `refkey10.txt`, `probes10.txt.gz` (~1.5 GB; request from author)

**Input:** Trimmed FASTQ files (`*_R1_tr.fastq.gz`, `*_R2_tr.fastq.gz`)  
**Output:** Per-sample species-level abundance CSV (read count + % RA)

**Resource note:** The k-mer classifier requires ~25 GB RAM.

---

### Step 4: Mitochondrial DNA Profiling

**What it does:** Identifies animal mitochondrial DNA in the metagenomes using two approaches: (a) the same k-mer classifier as Step 3 but with an mtDNA database (26,561 eukaryotic species), and (b) BURST v1.0 to quantify cattle-specific mtDNA reads (≥95% identity over 100 bp). Cattle mtDNA abundance is normalized to reads per million (RPM) for cross-sample comparison.

**Tools:**
- `nk10` with mtDNA database — full mtDNA profiling
- BURST v1.0 ([knights-lab/BURST](https://github.com/knights-lab/BURST)) — cattle mtDNA quantification at ≥95% identity over ≥100 bp

**Input:** Trimmed FASTQ files + mtDNA k-mer database + cattle mtDNA reference (for BURST)  
**Output:** mtDNA species profile per sample; cattle mtDNA RPM values

---

### Step 5: Microbial Source Tracking

**What it does:** Estimates the proportion of each air bacterial community derived from cattle feedlot pen soil using SourceTracker2, a Bayesian source-tracking tool. Source communities include 8 DREC pen soil samples collected in this study + 30 publicly available feedlot soil metagenomes from NCBI (BioProjects PRJNA292471 and PRJNA554493).

**Tools:**
- SourceTracker2 ([caporaso-lab/sourcetracker2](https://github.com/caporaso-lab/sourcetracker2)) — Bayesian source tracking
  - Parameters: 100 restarts, 10 draws per restart

**Input:** Species-level abundance table (from Step 3) for air samples (sinks) and pen soil samples (sources)  
**Output:** Per-sample source proportion estimates

---

### Step 6: Statistical & Comparative Analysis

**What it does:** Performs alpha/beta diversity analysis, differential abundance testing, and correlation analyses.

**Sub-steps:**
- **Alpha diversity:** Observed species richness (rarefied) and Shannon index computed at species level. Compared across gene-presence groups using Wilcoxon rank-sum test with Benjamini-Hochberg correction. (R v4.0.2, vegan v2.5-6)
- **Beta diversity:** Bray-Curtis dissimilarity matrix computed in Python 3.7.3. NMDS ordination via `sklearn.manifold`. Group significance tested with PERMANOVA (9,999 iterations) via `skbio`.
- **Differential abundance (LEfSe):** Identifies discriminatory taxa between virulence-gene-positive and -negative groups using default parameters.
- **Correlations:** Rank Biserial correlation for taxon RA vs. gene presence categories; Spearman rank for cattle mtDNA abundance vs. continuous variables.
- **Biomarker discovery (Pebblescout):** Used to query NCBI metagenomic datasets for *C. maris*-unique sequences to assess ecological specificity.

**Tools:**
- R v4.0.2 + vegan v2.5-6 — alpha diversity, rarefaction, Wilcoxon tests
- Python 3.7.3 + scikit-learn + scikit-bio — Bray-Curtis, NMDS, PERMANOVA
- LEfSe (Segata et al., 2011) — LDA effect size
- Pebblescout (Shiryev & Agarwala, 2024) — NCBI SRA k-mer search

---

## Snakemake Pipeline

### Configuration File (`config.yaml`)

```yaml
# === config.yaml ===
# Directories
raw_dir: "data/raw"
out_dir: "results"
db_dir: "databases"

# Samples (list sample prefixes here)
samples:
  - "air_sample_001"
  - "air_sample_002"
  - "pen_soil_001"
  # ... add all sample IDs

# Tool paths
trimmomatic_jar: "/path/to/trimmomatic-0.38.jar"
adapters: "/path/to/Trimmomatic-0.38/adapters/NexteraPE-PE.fa"
kmer_id_dir: "/path/to/kmer_id"
burst_bin: "/path/to/burst"

# Database paths
virulence_db: "databases/ecoli_virulence_serogroup_db"
cattle_mtdna_ref: "databases/cattle_mtdna.fasta"

# SourceTracker2
sourcetracker_sources: "metadata/sourcetracker_sources.txt"
sourcetracker_mapping: "metadata/sourcetracker_mapping.txt"

# Resource defaults
threads: 8
```

### Snakefile

```python
# === Snakefile ===
configfile: "config.yaml"

SAMPLES = config["samples"]
RAW = config["raw_dir"]
OUT = config["out_dir"]
DB  = config["db_dir"]
THR = config["threads"]

rule all:
    input:
        # Step 1: QC
        expand("{out}/qc/fastqc_raw/{s}_R1_fastqc.html", out=OUT, s=SAMPLES),
        expand("{out}/qc/fastqc_raw/{s}_R2_fastqc.html", out=OUT, s=SAMPLES),
        expand("{out}/trimmed/{s}_R1_tr.fastq.gz", out=OUT, s=SAMPLES),
        expand("{out}/trimmed/{s}_R2_tr.fastq.gz", out=OUT, s=SAMPLES),
        expand("{out}/qc/fastqc_trimmed/{s}_R1_tr_fastqc.html", out=OUT, s=SAMPLES),
        # Step 2: Virulence gene screening
        expand("{out}/virulence/{s}_blast_hits.txt", out=OUT, s=SAMPLES),
        # Step 3: Bacterial taxonomy
        expand("{out}/kmer_taxonomy/{s}_taxonomy.csv", out=OUT, s=SAMPLES),
        # Step 4: mtDNA profiling
        expand("{out}/mtdna/{s}_burst_cattle.txt", out=OUT, s=SAMPLES),
        expand("{out}/mtdna/{s}_mtdna_taxonomy.csv", out=OUT, s=SAMPLES),
        # Step 5: Source tracking
        f"{OUT}/sourcetracker/mixing_proportions.txt",
        # Step 6: Statistics
        f"{OUT}/stats/alpha_diversity.csv",
        f"{OUT}/stats/beta_diversity_nmds.csv",
        f"{OUT}/lefse/lefse_results.res",


# ============================================================
# STEP 1: READ QC & TRIMMING
# ============================================================

rule fastqc_raw:
    """Run FastQC on raw reads."""
    input:
        r1=f"{RAW}/{{sample}}_R1.fastq.gz",
        r2=f"{RAW}/{{sample}}_R2.fastq.gz",
    output:
        html1=f"{OUT}/qc/fastqc_raw/{{sample}}_R1_fastqc.html",
        html2=f"{OUT}/qc/fastqc_raw/{{sample}}_R2_fastqc.html",
    threads: 2
    resources:
        mem_mb=2000,
        time_min=30,
    shell:
        """
        mkdir -p {OUT}/qc/fastqc_raw
        fastqc -t {threads} -o {OUT}/qc/fastqc_raw {input.r1} {input.r2}
        """

rule trimmomatic:
    """Trim adapters (NexteraPE) and low-quality bases with Trimmomatic 0.38."""
    input:
        r1=f"{RAW}/{{sample}}_R1.fastq.gz",
        r2=f"{RAW}/{{sample}}_R2.fastq.gz",
    output:
        r1=f"{OUT}/trimmed/{{sample}}_R1_tr.fastq.gz",
        r1_unpaired=f"{OUT}/trimmed/{{sample}}_R1_unpaired.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2_tr.fastq.gz",
        r2_unpaired=f"{OUT}/trimmed/{{sample}}_R2_unpaired.fastq.gz",
    params:
        jar=config["trimmomatic_jar"],
        adapters=config["adapters"],
    threads: config["threads"]
    resources:
        mem_mb=8000,
        time_min=120,
    shell:
        """
        java -jar {params.jar} PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule fastqc_trimmed:
    """Run FastQC on trimmed reads for post-trim QC."""
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1_tr.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2_tr.fastq.gz",
    output:
        html1=f"{OUT}/qc/fastqc_trimmed/{{sample}}_R1_tr_fastqc.html",
    threads: 2
    resources:
        mem_mb=2000,
        time_min=30,
    shell:
        """
        mkdir -p {OUT}/qc/fastqc_trimmed
        fastqc -t {threads} -o {OUT}/qc/fastqc_trimmed {input.r1} {input.r2}
        """


# ============================================================
# STEP 2: E. COLI VIRULENCE GENE & SEROGROUP SCREENING
# ============================================================

rule make_blast_db:
    """Build BLAST database from custom E. coli virulence/serogroup gene sequences."""
    input:
        fasta=f"{DB}/ecoli_virulence_serogroup.fasta",
    output:
        ndb=f"{DB}/ecoli_virulence_serogroup_db.ndb",
    shell:
        """
        makeblastdb -in {input.fasta} \
            -dbtype nucl \
            -out {DB}/ecoli_virulence_serogroup_db
        """

rule virulence_blast:
    """
    Low-stringency BLAST of trimmed reads against virulence/serogroup database.
    Retain reads matching at >= 90% identity over >= 50 nt.
    This is the first-pass pre-filter; the paper then uses SeqMan NGen (commercial)
    for higher-stringency subtyping.
    """
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1_tr.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2_tr.fastq.gz",
        db=f"{DB}/ecoli_virulence_serogroup_db.ndb",
    output:
        hits=f"{OUT}/virulence/{{sample}}_blast_hits.txt",
    params:
        db_prefix=f"{DB}/ecoli_virulence_serogroup_db",
    threads: config["threads"]
    resources:
        mem_mb=16000,
        time_min=240,
    shell:
        """
        # Decompress and interleave reads for BLAST input
        mkdir -p {OUT}/virulence

        # Convert paired FASTQ to FASTA for BLAST
        zcat {input.r1} {input.r2} | \
            awk 'NR%4==1 {{printf ">%s\\n", substr($0,2)}} NR%4==2 {{print}}' \
            > {OUT}/virulence/{wildcards.sample}_reads.fasta

        # BLAST with low-stringency criteria: 90% identity, 50 nt alignment
        blastn -query {OUT}/virulence/{wildcards.sample}_reads.fasta \
            -db {params.db_prefix} \
            -out {output.hits} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -perc_identity 90 \
            -qcov_hsp_perc 33 \
            -num_threads {threads} \
            -evalue 1e-5

        # Filter for >= 50 nt alignment length
        awk -F'\\t' '$4 >= 50' {output.hits} > {output.hits}.tmp && \
            mv {output.hits}.tmp {output.hits}

        # Clean up intermediate FASTA
        rm -f {OUT}/virulence/{wildcards.sample}_reads.fasta
        """

# NOTE: The second-pass (SeqMan NGen v18) is commercial software.
# For an open-source alternative, you could map the pre-filtered reads
# back to the database with BWA-MEM or Bowtie2 at higher stringency.


# ============================================================
# STEP 3: BACTERIAL COMMUNITY PROFILING (k-mer taxonomy)
# ============================================================

rule compile_kmer_id:
    """Compile the kmer_id C++ classifier from the GitHub repo."""
    input:
        src=f"{config['kmer_id_dir']}/newkmer_10nx.cpp",
    output:
        bin=f"{config['kmer_id_dir']}/nk10",
    shell:
        """
        g++ -O3 {input.src} -o {output.bin} -lz
        """

rule kmer_taxonomy:
    """
    Run k-mer based bacterial taxonomy assignment.
    Requires bact10/ database subdirectory with bData10.txt, btree10.txt,
    refkey10.txt, probes10.txt.gz. Uses ~25 GB RAM.
    Input files must be named *_R1_tr.fastq.gz and *_R2_tr.fastq.gz.
    """
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1_tr.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2_tr.fastq.gz",
        nk10=f"{config['kmer_id_dir']}/nk10",
    output:
        csv=f"{OUT}/kmer_taxonomy/{{sample}}_taxonomy.csv",
    params:
        kmer_dir=config["kmer_id_dir"],
        input_dir=f"{OUT}/trimmed",
    threads: 1  # nk10 is single-threaded
    resources:
        mem_mb=28000,
        time_min=360,
    shell:
        """
        mkdir -p {OUT}/kmer_taxonomy

        # nk10 reads from a directory; it looks for *_R1_tr.fastq.gz / *_R2_tr.fastq.gz
        # Create a sample-specific symlink directory
        TMPDIR={OUT}/kmer_taxonomy/{wildcards.sample}_input
        mkdir -p $TMPDIR
        ln -sf $(realpath {input.r1}) $TMPDIR/
        ln -sf $(realpath {input.r2}) $TMPDIR/

        # Run kmer_id (requires bact10/ subdir to be in the kmer_id directory)
        cd {params.kmer_dir} && ./nk10 $TMPDIR/

        # Collect results with readbatch_10.py
        # (Adjust dir1 inside readbatch_10.py to point to results location)
        python3 {params.kmer_dir}/readbatch_10.py > {output.csv}

        rm -rf $TMPDIR
        """


# ============================================================
# STEP 4: MITOCHONDRIAL DNA PROFILING
# ============================================================

rule mtdna_kmer_profile:
    """
    Full mtDNA profiling using the kmer_id program with the mitochondrial
    database (26,561 eukaryotic species). Same tool as Step 3 but with
    mtDNA database files.
    """
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1_tr.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2_tr.fastq.gz",
        nk10=f"{config['kmer_id_dir']}/nk10",
    output:
        csv=f"{OUT}/mtdna/{{sample}}_mtdna_taxonomy.csv",
    params:
        kmer_dir=config["kmer_id_dir"],
    threads: 1
    resources:
        mem_mb=28000,
        time_min=360,
    shell:
        """
        mkdir -p {OUT}/mtdna

        # Similar to bacterial taxonomy, but ensure the mtDNA database
        # files are in the appropriate subdirectory within kmer_id.
        # The mitochondria_list.txt, mitochondria_data.txt, etc.
        # files in the repo define the mtDNA database structure.

        TMPDIR={OUT}/mtdna/{wildcards.sample}_input
        mkdir -p $TMPDIR
        ln -sf $(realpath {input.r1}) $TMPDIR/
        ln -sf $(realpath {input.r2}) $TMPDIR/

        cd {params.kmer_dir} && ./nk10 $TMPDIR/

        python3 {params.kmer_dir}/readbatch_10.py > {output.csv}
        rm -rf $TMPDIR
        """

rule burst_cattle_mtdna:
    """
    Quantify cattle mtDNA reads using BURST v1.0.
    Align trimmed reads to cattle mitochondrial genome reference.
    Filter for >= 95% identity over >= 100 bp.
    Compute reads per million (RPM) for normalization.
    """
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1_tr.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2_tr.fastq.gz",
    output:
        txt=f"{OUT}/mtdna/{{sample}}_burst_cattle.txt",
    params:
        burst=config["burst_bin"],
        ref=config["cattle_mtdna_ref"],
    threads: config["threads"]
    resources:
        mem_mb=16000,
        time_min=120,
    shell:
        """
        mkdir -p {OUT}/mtdna

        # Concatenate paired reads for BURST input
        zcat {input.r1} {input.r2} > {OUT}/mtdna/{wildcards.sample}_all.fastq

        # Run BURST with 95% identity threshold
        {params.burst} -r {params.ref} \
            -q {OUT}/mtdna/{wildcards.sample}_all.fastq \
            -o {OUT}/mtdna/{wildcards.sample}_burst_raw.txt \
            -i 0.95 \
            -t {threads}

        # Count total reads and cattle mtDNA hits (>= 100 bp alignment)
        TOTAL_READS=$(zcat {input.r1} | awk 'END{{print NR/4}}')
        TOTAL_READS=$((TOTAL_READS * 2))  # both pairs

        # Filter for alignment length >= 100 bp and compute RPM
        awk -v total="$TOTAL_READS" '
            $4 >= 100 {{count++}}
            END {{
                rpm = (count / total) * 1000000;
                print "sample\t{wildcards.sample}";
                print "total_reads\t" total;
                print "cattle_mtdna_reads\t" count;
                print "cattle_mtdna_rpm\t" rpm;
            }}
        ' {OUT}/mtdna/{wildcards.sample}_burst_raw.txt > {output.txt}

        rm -f {OUT}/mtdna/{wildcards.sample}_all.fastq
        """


# ============================================================
# STEP 5: SOURCE TRACKING (SourceTracker2)
# ============================================================

rule merge_abundance_tables:
    """
    Merge per-sample species abundance CSVs into a single OTU-style table
    for SourceTracker2 input.
    """
    input:
        csvs=expand(f"{OUT}/kmer_taxonomy/{{s}}_taxonomy.csv", s=SAMPLES),
    output:
        table=f"{OUT}/sourcetracker/abundance_table.txt",
    shell:
        """
        mkdir -p {OUT}/sourcetracker
        # Custom script to merge individual taxonomy CSVs into
        # a combined species x sample abundance matrix.
        python3 scripts/merge_taxonomy_tables.py \
            --input_dir {OUT}/kmer_taxonomy \
            --output {output.table}
        """

rule sourcetracker2:
    """
    Run SourceTracker2 with cattle pen soil samples as source
    and air samples as sink.
    Parameters: 100 restarts, 10 draws per restart.
    Source samples: 8 DREC + 30 NCBI (PRJNA292471 + PRJNA554493).
    """
    input:
        table=f"{OUT}/sourcetracker/abundance_table.txt",
        mapping=config["sourcetracker_mapping"],
    output:
        proportions=f"{OUT}/sourcetracker/mixing_proportions.txt",
    threads: config["threads"]
    resources:
        mem_mb=16000,
        time_min=480,
    shell:
        """
        sourcetracker2 gibbs \
            -i {input.table} \
            -m {input.mapping} \
            -o {OUT}/sourcetracker \
            --source_rarefaction_depth 0 \
            --sink_rarefaction_depth 0 \
            --restarts 100 \
            --draws_per_restart 10 \
            --jobs {threads}
        """


# ============================================================
# STEP 6: STATISTICAL ANALYSIS
# ============================================================

rule alpha_diversity:
    """
    Compute alpha diversity (observed species richness, Shannon index)
    at species level. Rarefied species richness computed in R with vegan.
    Wilcoxon rank-sum test with Benjamini-Hochberg correction.
    """
    input:
        table=f"{OUT}/sourcetracker/abundance_table.txt",
    output:
        csv=f"{OUT}/stats/alpha_diversity.csv",
    shell:
        """
        mkdir -p {OUT}/stats
        Rscript scripts/alpha_diversity.R \
            --input {input.table} \
            --output {output.csv}
        """

rule beta_diversity:
    """
    Compute Bray-Curtis dissimilarity and NMDS ordination.
    PERMANOVA with 9999 iterations for group comparisons.
    """
    input:
        table=f"{OUT}/sourcetracker/abundance_table.txt",
    output:
        nmds=f"{OUT}/stats/beta_diversity_nmds.csv",
    shell:
        """
        python3 scripts/beta_diversity.py \
            --input {input.table} \
            --output {output.nmds} \
            --permutations 9999
        """

rule lefse:
    """
    LEfSe analysis to identify discriminatory taxa between
    virulence-gene-positive and -negative groups.
    """
    input:
        table=f"{OUT}/sourcetracker/abundance_table.txt",
    output:
        res=f"{OUT}/lefse/lefse_results.res",
    shell:
        """
        mkdir -p {OUT}/lefse

        # Format input for LEfSe
        lefse_format_input.py \
            {input.table} \
            {OUT}/lefse/lefse_input.in \
            -c 1 -u 2 -o 1000000

        # Run LEfSe with default parameters
        lefse_run.py \
            {OUT}/lefse/lefse_input.in \
            {output.res}

        # Plot results
        lefse_plot_res.py \
            {output.res} \
            {OUT}/lefse/lefse_plot.png
        """
```

---

## Cluster Execution

### SLURM Profile (`profiles/slurm/config.yaml`)

```yaml
# Snakemake SLURM profile
executor: slurm
default-resources:
  slurm_partition: "high"
  mem_mb: 8000
  runtime: 120       # minutes
  cpus_per_task: 1
  slurm_account: "your_account"
jobs: 20
latency-wait: 60
```

### Bash Wrapper Script (`run_pipeline.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=leonard_meta
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err

# ── Load modules (adjust to your HPC) ──
module load fastqc/0.12.1
module load java/11           # for Trimmomatic
module load blast+/2.14.0
module load R/4.0.2
module load python/3.7

# ── Activate conda environment ──
source activate metagenomics  # or: conda activate metagenomics

# ── Create log directory ──
mkdir -p logs

# ── Run Snakemake with SLURM executor ──
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --profile profiles/slurm \
    --use-conda \
    --rerun-incomplete \
    --latency-wait 60 \
    --jobs 20 \
    -p
```

### Conda Environment (`envs/metagenomics.yaml`)

```yaml
name: metagenomics
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.7
  - fastqc=0.12.1
  - trimmomatic=0.38
  - blast=2.14.0
  - r-base=4.0.2
  - r-vegan=2.5.6
  - bioconductor-phyloseq
  - scikit-learn
  - scikit-bio
  - pandas
  - numpy
  - scipy
  - lefse
  - pip:
    - sourcetracker
```

---

## Notes on Reproducibility

1. **kmer_id database:** The bacterial k-mer probe database (`probes10.txt.gz`, ~1.5 GB) is too large for GitHub. Contact Mark Mammel (mmammel8) for access.

2. **SeqMan NGen v18** is commercial (DNASTAR). The BLAST pre-filter step is reproducible; for the high-stringency subtyping pass, consider BWA-MEM or Bowtie2 as open-source alternatives.

3. **BURST v1.0** is available at [github.com/knights-lab/BURST](https://github.com/knights-lab/BURST). It may require building from source.

4. **SourceTracker2** expects a QIIME-style mapping file with `SourceSink` column indicating whether each sample is a "source" (pen soil) or "sink" (air).

5. **Pebblescout** is an NCBI service for querying SRA datasets with k-mer probes; it runs through the NCBI web interface, not locally.

6. **Raw data** are deposited under NCBI BioProject PRJNA1252083.

7. **Correlation analyses** (Rank Biserial, Spearman) and the Pebblescout step are not included in the Snakemake workflow since they are downstream statistical analyses typically run interactively in R/Python notebooks.
