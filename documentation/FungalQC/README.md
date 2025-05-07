FungalQC.wdl
=========================

This WDL workflow performs de novo genome assembly, taxonomic identification, and quality control of paired-end eukaryotic NGS data.\
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


Workflow Outputs
-------

### Version Information

| Output | Type | Description |
| --- | --- | --- |
| theiaeuk_illumina_pe_version | String | Pipeline version |
| theiaeuk_illumina_pe_analysis_date | String | Analysis date |

### Read Subsampling (RASUSA)

| Output | Type | Description |
| --- | --- | --- |
| rasusa_version | String | RASUSA version |
| read1_subsampled | File | Subsampled forward reads |
| read2_subsampled | File | Subsampled reverse reads |

### Read Metadata

| Output | Type | Description |
| --- | --- | --- |
| seq_platform | String | Sequencing platform |

### Read Screening

| Output | Type | Description |
| --- | --- | --- |
| read_screen_raw | String | Raw read screening result |
| read_screen_raw_tsv | File | Raw read screening details |
| read_screen_clean | String | Clean read screening result |
| read_screen_clean_tsv | File | Clean read screening details |

### Read QC (fastq_scan)

| Output | Type | Description |
| --- | --- | --- |
| fastq_scan_num_reads_raw1 | Int | Number of raw forward reads |
| fastq_scan_num_reads_raw2 | Int | Number of raw reverse reads |
| fastq_scan_num_reads_raw_pairs | String | Number of raw read pairs |
| fastq_scan_version | String | fastq_scan version |
| fastq_scan_num_reads_clean1 | Int | Number of clean forward reads |
| fastq_scan_num_reads_clean2 | Int | Number of clean reverse reads |
| fastq_scan_num_reads_clean_pairs | String | Number of clean read pairs |
| fastq_scan_raw1_json | File | Raw forward reads JSON report |
| fastq_scan_raw2_json | File | Raw reverse reads JSON report |
| fastq_scan_clean1_json | File | Clean forward reads JSON report |
| fastq_scan_clean2_json | File | Clean reverse reads JSON report |

### Read Trimming (Trimmomatic)

| Output | Type | Description |
| --- | --- | --- |
| trimmomatic_version | String | Trimmomatic version |
| trimmomatic_docker | String | Trimmomatic Docker image |

### Read QC (FastQC)

| Output | Type | Description |
| --- | --- | --- |
| fastqc_num_reads_raw1 | Int | Number of raw forward reads |
| fastqc_num_reads_raw2 | Int | Number of raw reverse reads |
| fastqc_num_reads_raw_pairs | String | Number of raw read pairs |
| fastqc_num_reads_clean1 | Int | Number of clean forward reads |
| fastqc_num_reads_clean2 | Int | Number of clean reverse reads |
| fastqc_num_reads_clean_pairs | String | Number of clean read pairs |
| fastqc_raw1_html | File | Raw forward reads HTML report |
| fastqc_raw2_html | File | Raw reverse reads HTML report |
| fastqc_clean1_html | File | Clean forward reads HTML report |
| fastqc_clean2_html | File | Clean reverse reads HTML report |
| fastqc_version | String | FastQC version |
| fastqc_docker | String | FastQC Docker image |

### Read QC (fastp)

| Output | Type | Description |
| --- | --- | --- |
| fastp_version | String | fastp version |
| fastp_html_report | File | fastp HTML report |

### Cleaned Reads

| Output | Type | Description |
| --- | --- | --- |
| bbduk_docker | String | BBDuk Docker image |
| read1_clean | File | Cleaned forward reads |
| read2_clean | File | Cleaned reverse reads |

### Read QC (CG Pipeline)

| Output | Type | Description |
| --- | --- | --- |
| r1_mean_q_raw | Float | Raw forward reads mean quality |
| r2_mean_q_raw | Float | Raw reverse reads mean quality |
| combined_mean_q_raw | Float | Raw combined reads mean quality |
| combined_mean_q_clean | Float | Clean combined reads mean quality |
| r1_mean_readlength_raw | Float | Raw forward reads mean length |
| r2_mean_readlength_raw | Float | Raw reverse reads mean length |
| combined_mean_readlength_raw | Float | Raw combined reads mean length |
| combined_mean_readlength_clean | Float | Clean combined reads mean length |

### Assembly (Shovill)

| Output | Type | Description |
| --- | --- | --- |
| assembly_fasta | File | Assembled contigs |
| contigs_gfa | File | Assembly graph in GFA format |
| contigs_fastg | File | Assembly graph in FASTG format |
| contigs_lastgraph | File | Assembly graph in LastGraph format |
| shovill_pe_version | String | Shovill version |

### EukCC

| Output | Type | Description |
| --- | --- | --- |
| EukCC_output_folder | File | EukCC output directory |

### Assembly QC (QUAST)

| Output | Type | Description |
| --- | --- | --- |
| quast_report | File | QUAST report |
| quast_version | String | QUAST version |
| assembly_length | Int | Total assembly length |
| number_contigs | Int | Number of contigs |
| n50_value | Int | N50 value |
| quast_gc_percent | Float | GC content percentage |

### Assembly QC (CG Pipeline)

| Output | Type | Description |
| --- | --- | --- |
| cg_pipeline_report_raw | File | CG Pipeline report for raw reads |
| cg_pipeline_docker | String | CG Pipeline Docker image |
| est_coverage_raw | Float | Estimated coverage from raw reads |
| cg_pipeline_report_clean | File | CG Pipeline report for clean reads |
| est_coverage_clean | Float | Estimated coverage from clean reads |

### Assembly QC (BUSCO)

| Output | Type | Description |
| --- | --- | --- |
| busco_version | String | BUSCO version |
| busco_docker | String | BUSCO Docker image |
| busco_database | String | BUSCO database used |
| busco_results | String | BUSCO results summary |
| busco_report | File | BUSCO report |

### Taxonomic ID (GAMBIT)

| Output | Type | Description |
| --- | --- | --- |
| gambit_report | File | GAMBIT report |
| gambit_closest_genomes | File | GAMBIT closest genomes |
| gambit_predicted_taxon | String | Predicted taxonomic classification |
| gambit_predicted_taxon_rank | String | Taxonomic rank of prediction |
| gambit_version | String | GAMBIT version |
| gambit_db_version | String | GAMBIT database version |
| gambit_docker | String | GAMBIT Docker image |

### Taxonomic ID (Kraken2)

| Output | Type | Description |
| --- | --- | --- |
| kraken2_report | File | Kraken2 report |

### QC Check

| Output | Type | Description |
| --- | --- | --- |
| qc_check | String | QC check result |
| qc_standard | File | QC standards used |

### Species-Specific Analysis

| Output | Type | Description |
| --- | --- | --- |
| cladetyper_clade | String | Identified clade |
| cladetyper_gambit_version | String | CladeTyper GAMBIT version |
| cladetyper_docker_image | String | CladeTyper Docker image |
| cladetyper_annotated_reference | String | Annotated reference used |
| theiaeuk_snippy_variants_query | String | Snippy variants query |
| theiaeuk_snippy_variants_query_check | String | Snippy variants query check |
| theiaeuk_snippy_variants_hits | String | Snippy variants hits |
| theiaeuk_snippy_variants_gene_query_results | String | Snippy variants gene query results |

## Task Descriptions and Enhancements

### Task: [gambit](https://github.com/broadinstitute/idmp-fungal-pipelines/blob/main/tasks/taxon_id/task_gambit.wdl)

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
- `gambit_predicted_taxon_rank`: Rank of the predicted taxon
- `gambit_next_taxon`: Next closest species match
- `gambit_next_taxon_rank`: Rank of the next closest match
- `gambit_version`: GAMBIT software version
- `gambit_db_version`: Version of the GAMBIT database used
- `merlin_tag`: Filter tag used to control downstream steps
- `gambit_docker`: Docker image used to run the GAMBIT task

**Enhancements:**
A custom Python block was added to enforce species-specific gating:

```python
# Write out warning if the predicted taxon is not the expected taxon
gambit_expected_taxon = "~{gambit_expected_taxon}"
if merlin_tag != gambit_expected_taxon:
    print(f"WARNING! Pipeline is configured to only proceed for {gambit_expected_taxon}. Found: {merlin_tag}", file=sys.stderr)
```

This allows the overall FungalQC workflow to end early if the sample is not identified as Candidozyma auris, preventing unnecessary downstream processing for non-target species.

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

Additional Notes
----------------

-   The workflow is designed for eukaryotic genomes, with default parameters optimized for fungal genomes.
-   The default expected taxon is set to "Candidozyma auris" (Candida auris).
-   The workflow includes validation steps that can be turned off with the `skip_screen` parameter.
-   Species-specific analysis is only performed when the identified organism matches the expected taxon.

References
----------

-   [TheiaEuk Pipeline Documentation](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaeuk/)

-   [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)

-   [EukCC GitHub](https://github.com/Finn-Lab/EukCC)