FungalQC.wdl
=========================

This WDL workflow performs de novo genome assembly, taxonomic identification, and quality control of paired-end eukaryotic NGS data.\
It is based on the [Theiagen TheiaEuk pipeline](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/), which provides detailed descriptions of the standard tasks and processes used.

Workflow Overview
-----------------

The workflow follows the overall structure of the TheiaEuk pipeline, including steps for read QC, assembly, taxonomic identification, and assembly quality control.

For full descriptions of the standard tasks, refer to [Theiagen's documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/).

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

### `bam_filter_fixmates` Task (version 1.0)

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

References
----------

-   [TheiaEuk Pipeline Documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/)

-   [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

-   [EukCC GitHub](https://github.com/Finn-Lab/EukCC)