FungalQC.wdl
=========================

This WDL workflow performs de novo genome assembly, taxonomic identification, and quality control of paired-end eukaryotic NGS data.

It is based on the [Theiagen TheiaEuk pipeline](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/), which provides detailed descriptions of the standard tasks and processes used.

Workflow Overview
-----------------

The workflow follows the overall structure of the TheiaEuk pipeline, including steps for read QC, assembly, taxonomic identification, and assembly quality control.
For full descriptions of the standard tasks, refer to [Theiagen's documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/).

The FungalQc.wdl combines multiple tools to perform:

1.  Read quality control and trimming
2.  Subsampling for high-coverage datasets
3.  De-novo genome assembly
4.  Assembly quality assessment
5.  Taxonomic identification
6.  Fungal species-specific analysis

Workflow Inputs
------

### Required Inputs

| Input | Type | Description |
| --- | --- | --- |
| samplename | String | Name for the sample |
| read1 | File | Forward reads FASTQ file |
| read2 | File | Reverse reads FASTQ file |

### Optional Configuration Inputs

| Input | Type | Description | Default |
| --- | --- | --- | --- |
| seq_method | String | Sequencing method | "ILLUMINA" |
| call_rasusa | Boolean | Whether to subsample reads | true |
| min_reads | Int | Minimum number of reads required | 30000 |
| min_basepairs | Int | Minimum number of base pairs required | 45000000 |
| min_genome_length | Int | Minimum expected genome length | 9000000 |
| max_genome_length | Int | Maximum expected genome length | 178000000 |
| min_coverage | Int | Minimum required coverage | 10 |
| min_proportion | Int | Minimum proportion required | 40 |
| trim_min_length | Int | Minimum read length after trimming | 75 |
| trim_quality_min_score | Int | Minimum quality score for trimming | 20 |
| trim_window_size | Int | Window size for trimming | 10 |
| busco_memory | Int | Memory allocation for BUSCO (GB) | 24 |
| busco_docker_image | String | Docker image for BUSCO | "us-docker.pkg.dev/general-theiagen/ezlabgva/busco.3.2_cv1" |
| skip_screen | Boolean | Skip read screening step | false |
| qc_check_table | File? | Optional QC check table |  |
| expected_taxon | String? | Expected taxonomic classification |  |
| genome_length | Int? | Expected genome length |  |
| subsample_coverage | Float | Target coverage for subsampling | 150 |
| cpu | Int | Number of CPUs to use | 8 |
| memory | Int | Memory allocation (GB) | 16 |
| gambit_db_genomes | File | GAMBIT database genomes | gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-metadata-1.0.0-20241213.gdb |
| gambit_db_signatures | File | GAMBIT database signatures | gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-signatures-1.0.0-20241213.gs |
| contamination_percent_threshold | Float | EukCC contamination threshold | 5.0 |
| gambit_expected_taxon | String | Expected taxonomic classification for GAMBIT | "Candidozyma auris" |

---

### Task: [kraken2](https://github.com/broadinstitute/idmp-fungal-pipelines/blob/main/tasks/taxon_id/task_kraken2.wdl)

**Purpose:**
Adds an additional taxonomic identification method alongside GAMBIT.

**Inputs:**
- `assembly`: Contig file for classification
- `kraken2_db`: Optional Path to Kraken2 database. If left empty, Kraken2 will use the database baked into the [Docker image](https://github.com/broadinstitute/idmp-fungal-pipelines/blob/main/dockers/kraken2/Dockerfile).
- `read1`: Paired-end read file 1
- `read2`:  Paired-end read file 2

**Outputs:**
- `kraken2_report`: Classification report
- `classified_reads`: File with classified reads
- `kraken2_report_taxon_name`: Predicted taxon name (species level)


---

### Task: [EukCC](https://github.com/broadinstitute/idmp-fungal-pipelines/blob/main/tasks/quality_control/advanced_metrics/task_EukCC.wdl)

**Purpose:**
Assess genome completeness and contamination for eukaryotic assemblies.

**Inputs:**
- `assembly`: Assembled genome file
- `eukcc_db_path`: Optional Path to EukCC database. If left empty, EukCC will use the database baked into the [Docker image](https://github.com/broadinstitute/idmp-fungal-pipelines/blob/main/dockers/eukcc/Dockerfile).
- `contamination_percent_threshold`: Threshold for allowable contamination

**Outputs:**
- `eukcc_csv`: Full EukCC result CSV file
- `completeness`: Parsed completeness percentage
- `contamination`: Parsed contamination percentage

---

### Task: [bam_filter_fixmates](https://github.com/broadinstitute/idmp-fungal-pipelines/blob/main/tasks/quality_control/read_filtering/task_filter_bam.wdl)

**Purpose:**
This task filters a BAM file to retain only properly paired reads and indexes the resulting BAM for downstream processing.

**Inputs:**
- `input_bam`: Input BAM file

**Outputs:**
- `filtered_bam`:  Filtered BAM file
- `filtered_bai`: Filtered BAM index file

**Role in Workflow:**
This task is part of `theiaeuk_merlin_typing.wdl`, which is a workflow focused on *Candida auris* typing. Filtering for properly paired reads ensures accurate variant calling and typing within the *C. auris* analysis pipeline by removing ambiguous or low-quality alignments.

---


### Merlin Typing WDL Refinement

-   **Renamed File:** `merlin_magic.wdl` has been renamed to `theiaeuk_merlin_typing.wdl`.

-   **Content Simplification:** All functionality unrelated to *Candida auris* processing has been removed. The resulting WDL is streamlined to focus exclusively on *C. auris* species detection and typing workflows.


**Enhancements:**
A custom Python block was added to enforce species-specific gating:

```python
# Write out warning if the predicted taxon is not the expected taxon
gambit_expected_taxon = "~{gambit_expected_taxon}"
if merlin_tag != gambit_expected_taxon:
    print(f"WARNING! Pipeline is configured to only proceed for {gambit_expected_taxon}. Found: {merlin_tag}", file=sys.stderr)
```

This allows the overall FungalQC workflow to end early if the sample is not identified as Candidozyma auris, preventing unnecessary downstream processing for non-target species.


Additional Notes
----------------

-   The workflow is designed for eukaryotic genomes, with default parameters optimized for fungal genomes.
-   The default expected taxon is set to "Candidozyma auris" (Candida auris).
-   Species-specific analysis is only performed when the identified organism matches the expected taxon.

References
----------

-   [TheiaEuk Pipeline Documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/)

-   [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

-   [EukCC GitHub](https://github.com/Finn-Lab/EukCC)