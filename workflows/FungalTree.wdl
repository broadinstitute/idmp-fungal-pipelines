import "../workflows/FungalVariantCallingGatk3.wdl" as FungalVariantCallingGatk3

workflow FungalTree {

    meta {
        description: "FungalTree is a WDL-based pipeline for variant calling and phylogenetic analysis in fungal haploid genomes, using a reference GenBank file and aligned BAM files to generate filtered variant calls and a maximum-likelihood phylogenetic tree via IQ-TREE2."
        allowNestedInputs: "true"
    }

    ## config params
    # input data
    String analysis_name
    File ref_gbff
    Array[String] input_samples
    Array[File] input_bams

    # hard filtering params: both of these params are required
    String snp_filter_expr
    String indel_filter_expr

    # IQ-TREE2 parameters
    String? iqtree2_model
    Int iqtree2_bootstraps = 1000
    Int alrt = 1000
    String? iqtree2_opts

    call GbffToFasta {
        input:
            ref_gbff = ref_gbff
    }

    call GenerateRefFiles {
        input:
            ref_fasta = GbffToFasta.reference_fasta,
            input_bam = input_bams[0]
    }

    call FungalVariantCallingGatk3.FungalVariantCallingGatk3 as FungalVariantCallingGatk3 {
    input:
           analysis_name = analysis_name,
           input_samples = input_samples,
           input_bams = input_bams,
           reference_fasta = GenerateRefFiles.reference_fasta,
           ref_dict = GenerateRefFiles.ref_dict,
           ref_index = GenerateRefFiles.ref_index,
           snp_filter_expr = snp_filter_expr ,
           indel_filter_expr = indel_filter_expr
    }


    call VCFToFasta {
        input:
        vcf_file = FungalVariantCallingGatk3.hard_filtered_gvcf
    }

    call IqTree2 {
        input:
        alignment = VCFToFasta.alignment_fasta,
        iqtree2_model = iqtree2_model,
        iqtree2_bootstraps = iqtree2_bootstraps,
        alrt = alrt,
        iqtree2_opts = iqtree2_opts,
        cluster_name = analysis_name
    }

    output {
        File hard_filtered_gvcf = FungalVariantCallingGatk3.hard_filtered_gvcf
        File IqTree2_ml_tree = IqTree2.ml_tree
        String IqTree2_run_date = IqTree2.date
        String IqTree2_model_used = IqTree2.iqtree2_model_used
    }
}

task GbffToFasta {
    File ref_gbff
    String ref_gbff_basename = basename(ref_gbff, ".fasta")

    Int disk_size = 50
    Int mem_size_gb = 16
    String docker = "python:3.11-slim"

    command <<<
        set -euo pipefail
        pip install biopython
        python <<CODE

        from Bio import SeqIO
        input_file = "${ref_gbff}"
        output_file = "${ref_gbff_basename}.fasta"
        with open(output_file, "w") as out_fasta:
            for record in SeqIO.parse(input_file, "genbank"):
                 SeqIO.write(record, out_fasta, "fasta")
        CODE

        # strip the .1 in the reference sequence name
        sed -i -E 's/^(>[^ ]+)\.[0-9]+/\1/' "${ref_gbff_basename}.fasta"
    >>>
    output {
        File reference_fasta =  "${ref_gbff_basename}.fasta"
    }
    runtime {
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"

    }
}

task GenerateRefFiles {
    File ref_fasta
    File input_bam
    String ref_fasta_basename = basename(ref_fasta, ".fasta")

    Int disk_size = 50
    Int mem_size_gb = 16
    String docker = "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.0-0.7.15-2.26.3-1634165082"

    command <<<
        set -euo pipefail

        echo ${ref_fasta_basename}

        cp ${ref_fasta} ${ref_fasta_basename}.fasta
        /usr/gitc/bwa index ${ref_fasta_basename}.fasta

        java -Xms1000m -Xmx1000m  -jar /usr/gitc/picard.jar CreateSequenceDictionary R=${ref_fasta_basename}.fasta O=${ref_fasta_basename}.dict

        samtools faidx ${ref_fasta_basename}.fasta


    >>>
    output {
    File ref_sa = "${ref_fasta_basename}.fasta.sa"
    File ref_bwt = "${ref_fasta_basename}.fasta.bwt"
    File ref_amb = "${ref_fasta_basename}.fasta.amb"
    File ref_ann = "${ref_fasta_basename}.fasta.ann"
    File ref_pac = "${ref_fasta_basename}.fasta.pac"
    File ref_dict = "${ref_fasta_basename}.dict"
    File ref_index = "${ref_fasta_basename}.fasta.fai"
    File reference_fasta = "${ref_fasta_basename}.fasta"

    }
    runtime {
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"

    }
}

task IqTree2 {
    File alignment
    String cluster_name
    String? iqtree2_model # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
    Int iqtree2_bootstraps = 1000 #  Ultrafast bootstrap replicates
    Int alrt = 1000 # SH-like approximate likelihood ratio test (SH-aLRT) replicates
    String? iqtree2_opts

    String docker = "us-docker.pkg.dev/general-theiagen/staphb/iqtree2:2.1.2"
    Int disk_size = 50
    Int cpu = 4
    Int memory = 32

  command <<<
        # Date and version control
        date | tee DATE

        # Get IQ-TREE version
        iqtree2 --version | grep version | sed 's|.*version|version|;s| COVID-edition for Linux.*||' | tee VERSION

        # Check if iqtree2_model input is set and output for sanity
        if [ -n "${iqtree2_model}" ]; then
            echo "DEBUG: User provided iqtree2_model ${iqtree2_model}, will use this for running iqtree2"
            IQTREE2_MODEL="${iqtree2_model}"
        else
            echo "DEBUG: User did not supply an iqtree2_model input, will use iqtree2's model finder"
        fi

        # Sanity check
        echo "DEBUG: IQTREE2_MODEL is set to: " $IQTREE2_MODEL

        # Make sure there are more than 3 genomes in the dataset
        numGenomes=$(grep -o '>' ${alignment} | wc -l)
        if [ "$numGenomes" -gt 3 ]; then
            cp ${alignment} ./msa.fasta

            # Run iqtree2
            if [[ -v IQTREE2_MODEL ]] ; then # iqtree2 model set; use -m tag
                echo "DEBUG: running iqtree2 with the -m flag which is used to provide a model; user-specified " $IQTREE2_MODEL
                iqtree2 \
                    -nt AUTO \
                    -s msa.fasta \
                    -m $IQTREE2_MODEL \
                    -bb ${iqtree2_bootstraps} \
                    -alrt ${alrt} ${iqtree2_opts}

                # Write the iqtree2_model used to a txt file for output as a string
                echo $IQTREE2_MODEL | tee IQTREE2_MODEL.TXT
            else # iqtree model is not set; do not use -m tag
                echo "DEBUG: running iqtree2 without the -m flag which is used to provide a model. Will default to iqtree2 default (Model Finder)"
                iqtree2 \
                    -nt AUTO \
                    -s msa.fasta \
                    -bb ${iqtree2_bootstraps} \
                    -alrt ${alrt} ${iqtree2_opts}

                # Determine iqtree2_model used by parsing log file
                grep "Best-fit model" msa.fasta.log | sed 's|Best-fit model: ||g;s|chosen.*||' | tee IQTREE2_MODEL.TXT
            fi

            # Rename the final output newick file
            cp -v msa.fasta.contree ${cluster_name}_iqtree.nwk
        else
            echo "ERROR: Not enough genomes provided; more than 3 are required to run iqtree2"
            exit 1
        fi
    >>>
  output {
    String date = read_string("DATE")
    String iqtree2_version = read_string("VERSION")
    File ml_tree = "${cluster_name}_iqtree.nwk"
    String iqtree2_model_used = read_string("IQTREE2_MODEL.TXT")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}

task VCFToFasta {
    File vcf_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")

    Int cpu = 4
    Int disk_size = 50
    Int memory = 16
  command <<<
    python3 /app/vcf2matrix.py \
      -f \
      -i ${vcf_file} \
      --output-prefix ${vcf_basename}
  >>>
  output {
    File alignment_fasta = "${vcf_basename}.min4.fasta"
  }
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/vcftomsa:1.0.0"
        memory: memory + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB"
        preemptible: 0
    }
}
