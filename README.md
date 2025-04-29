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

Requirements
------------

-   Paired-end Illumina sequencing data (FASTQ format)

-   Reference databases for Kraken2 and GAMBIT

-   A WDL execution engine such as Cromwell

Outputs
-------

The workflow generates:

-   Assembled genome FASTAs

-   Taxonomic reports from GAMBIT and Kraken2

-   Assembly quality metrics from EukCC, BUSCO, and QUAST

-   *C. auris*-specific typing results (if applicable)

References
----------

-   [TheiaEuk Pipeline Documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/)

-   [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

-   [EukCC GitHub](https://github.com/Finn-Lab/EukCC)