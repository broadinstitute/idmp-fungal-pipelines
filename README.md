FungalQC.wdl
=========================

This WDL workflow performs de novo genome assembly, taxonomic identification, and quality control of paired-end eukaryotic NGS data.\
It is based on the [Theiagen TheiaEuk pipeline](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/), which provides detailed descriptions of the standard tasks and processes used.

Workflow Overview
-----------------

The workflow follows the overall structure of the TheiaEuk pipeline, including steps for read QC, assembly, taxonomic identification, and assembly quality control.

For full descriptions of the standard tasks, refer to [Theiagen's documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/).


## Task Descriptions and Enhancements

### Task: gambit

**Purpose:**
Taxonomic identification of assembled genome using GAMBIT.

**Inputs:**
- `assembly`: Assembly file to classify
- `samplename`: Sample ID used in naming outputs
- `gambit_db_genomes` and `gambit_db_signatures`: Optional custom GAMBIT database files
- `gambit_expected_taxon`: Target taxon (e.g., "Candidozyma auris") for gating

**Outputs:**
- `gambit_report_file`: Full JSON GAMBIT report
- `gambit_closest_genomes_file`: CSV of closest genome matches
- `gambit_predicted_taxon`: Predicted species name
- `merlin_tag`: Filter tag used to control downstream steps

**Enhancements:**
A custom Python block was added to enforce species-specific gating:

```python
# Write out warning if the predicted taxon is not the expected taxon
gambit_expected_taxon = "~{gambit_expected_taxon}"
if merlin_tag != gambit_expected_taxon:
    print(f"WARNING! Pipeline is configured to only proceed for {gambit_expected_taxon}. Found: {merlin_tag}", file=sys.stderr)
```

This allows the overall FungalQC workflow to fail early if the sample is not identified as Candidozyma auris, preventing unnecessary downstream processing for non-target species.

---

### Task: kraken2

**Purpose:**
Adds an additional taxonomic identification method alongside GAMBIT.

**Inputs:**
- `assembly`: Contig file for classification
- `kraken2_db`: Optional Path to Kraken2 database
- `read1`:
- `read2`:

**Outputs:**
- `kraken2_report`: Classification report
- `classified_reads`:
- `kraken2_report_taxon_name`:


---

### Task: EukCC

**Purpose:**
Assess genome completeness and contamination for eukaryotic assemblies.

**Inputs:**
- `assembly`: Assembled genome file
- `eukcc_db_path`: Optional Path to EukCC database
- `contamination_percent_threshold`:

**Outputs:**
- `eukcc_csv`:
- `completeness`:
- `contamination`:

---

### Task: quast

**Purpose:**
Evaluate assembly metrics including N50, total contig length, and number of contigs.

**Inputs:**
- `assembly`: Assembled genome file

**Outputs:**
- `quast_report`: Summary of assembly statistics

**Enhancements:**
- No enhancements or parameter changes. Used for basic QC.

---







Modifications and Enhancements
------------------------------

This version of the workflow introduces several key changes and additions:

### EukCC Integration

-   **Purpose:** Enhances genome quality evaluation for eukaryotic assemblies.

-   **Description:** EukCC assesses assembly completeness and contamination using lineage-specific markers, offering more accurate and detailed QC metrics specifically tuned for eukaryotic organisms.

### Kraken2 Integration

-   **Purpose:** Adds an additional taxonomic identification method alongside GAMBIT.

-   **Description:** Kraken2 is a high-speed k-mer-based classifier that enables broader detection of sample contamination and unexpected species. It complements GAMBIT's targeted taxonomic assignment by providing full-sample classification.

### Merlin Typing WDL Refinement

-   **Renamed File:** `merlin_magic.wdl` has been renamed to `theiaeuk_merlin_typing.wdl`.

-   **Content Simplification:** All functionality unrelated to *Candida auris* processing has been removed. The resulting WDL is streamlined to focus exclusively on *C. auris* species detection and typing workflows.

Task Descriptions
-----------------

### `bam_filter_fixmates` Task

**Purpose:**
This task filters a BAM file to retain only properly paired reads and indexes the resulting BAM for downstream processing.

**Description:**
`bam_filter_fixmates` takes a BAM file as input and produces a filtered BAM containing only reads with the SAM flag `0x2` (properly paired). The task also creates indexes for both the original and filtered BAM files.

**Inputs:**

- `input_bam`: Input BAM file
- `output_prefix`: Prefix for output files
- `cpu`: Number of CPUs to use (default: 2)
- `memory`: Memory in GB (default: 8)
- `docker`: Docker image for `samtools` (default: `staphb/samtools:1.19`)

**Outputs:**

- Filtered BAM file (`*_filtered.bam`)
- BAM index file (`*_filtered.bam.bai`)

**Role in Workflow:**
This task is part of `theiaeuk_merlin_typing.wdl`, which is a refined module focused on *Candida auris* typing. Filtering for properly paired reads ensures accurate variant calling and typing within the *C. auris* analysis pipeline by removing ambiguous or low-quality alignments.

### `gambit` Task

**Purpose:**
Performs taxonomic identification using GAMBIT on an assembled genome and produces detailed species-level classification, along with closest genome matches.

**Inputs:**

- `assembly`: Input genome assembly FASTA
- `samplename`: Sample identifier
- `docker`: GAMBIT Docker image
- `gambit_db_genomes`, `gambit_db_signatures`: GAMBIT reference databases
- `gambit_expected_taxon`: Expected taxon name, e.g., `"Candidozyma auris"`

**Outputs:**

- JSON report (`*_gambit.json`)
- CSV of closest genome matches (`*_gambit_closest.csv`)
- Predicted taxon information and version metadata
- `merlin_tag`: Canonical tag used to gate further workflow execution

**Enhancement:**
A custom Python block was added to enforce species-specific gating:
```python
# Write out warning if the predicted taxon is not the expected taxon
gambit_expected_taxon = "~{gambit_expected_taxon}"
if merlin_tag != gambit_expected_taxon:
    print(f"WARNING! Pipeline is configured to only proceed for {gambit_expected_taxon}.
```

This allows the overall FungalQC workflow to fail early if the sample is not identified as Candidozyma auris, preventing unnecessary downstream processing for non-target species.

References
----------

-   [TheiaEuk Pipeline Documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/)

-   [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

-   [EukCC GitHub](https://github.com/Finn-Lab/EukCC)