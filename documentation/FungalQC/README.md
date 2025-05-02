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

### Task: bam_filter_fixmates

**Purpose:**
This task filters a BAM file to retain only properly paired reads and indexes the resulting BAM for downstream processing.

**Inputs:**
- `input_bam`: Input BAM file

**Outputs:**
- `filtered_bam`:  Filtered BAM file
- `filtered_bai`: Filtered BAM index file

**Role in Workflow:**
This task is part of `theiaeuk_merlin_typing.wdl`, which is a refined module focused on *Candida auris* typing. Filtering for properly paired reads ensures accurate variant calling and typing within the *C. auris* analysis pipeline by removing ambiguous or low-quality alignments.

---


### Merlin Typing WDL Refinement

-   **Renamed File:** `merlin_magic.wdl` has been renamed to `theiaeuk_merlin_typing.wdl`.

-   **Content Simplification:** All functionality unrelated to *Candida auris* processing has been removed. The resulting WDL is streamlined to focus exclusively on *C. auris* species detection and typing workflows.


References
----------

-   [TheiaEuk Pipeline Documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/)

-   [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

-   [EukCC GitHub](https://github.com/Finn-Lab/EukCC)