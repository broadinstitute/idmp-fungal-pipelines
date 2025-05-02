## Copyright Broad Institute 2020
##
## Variant Calling pipeline for fungal haploid genomes
## Developed by Xiao Li (xiaoli@broadinstitute.org)
## Fungal Genomics Group, Infectious Disease and Microbiome Program.
## The Broad Institute of MIT and Harvard
##
## Verions of the software used in this WDL:
##   PICARD_VER=1.782
##   GATK37_VER=3.7-93-ge9d8068
##   SAMTOOLS_VER=1.3.1
##   BWA_VER=0.7.12
##   TABIX_VER=0.2.5_r1005
##   BGZIP_VER=1.3
##
## Cromwell version support
## - Successfully tested on v28, v30 and v36
## - Does not work on versions < v23 due to output syntax
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.


workflow FungalTree {
    ## config params
    # input data
    String run_name

    File ref
    File ref_gbff
    #File ref_sa
    #File ref_bwt
    #File ref_amb
    #File ref_ann
    #File ref_pac
    #File ref_dict
    #File ref_index
    Array[String] input_samples
    Array[File] input_bams

    # mem size/ disk size params
    Int small_mem_size_gb
    Int med_mem_size_gb
    Int large_mem_size_gb
    Int extra_large_mem_size_gb

    Int disk_size
    Int med_disk_size
    Int large_disk_size
    Int extra_large_disk_size

    String docker

    String picard_path
    String gatk_path

    # whether perform alignment or not to the input BAM files
    Boolean do_align

    # hard filtering params: both of these params are required
    String snp_filter_expr
    String indel_filter_expr

    # IQ-TREE2 parameters
    String? iqtree2_model
    Int iqtree2_bootstraps = 1000
    Int alrt = 1000
    String? iqtree2_opts


    ## task calls
    call GbffToFasta {
        input:
            ref_gbff = ref_gbff
    }

    call GenerateRefFiles {
        input:
            ref_fasta = GbffToFasta.reference_fasta,
            input_bam = input_bams[0] # Just need one BAM for this?
    }

    # run pipeline on each sample, in parallel
    scatter(i in range(length(input_samples))) {
        String sample_name = input_samples[i]
        String input_bam = input_bams[i]

        if (do_align) {
            call SamToFastq {
                input:
                in_bam = input_bam,
                sample_name = sample_name,
                disk_size = large_disk_size,
                mem_size_gb = small_mem_size_gb,
                docker = docker,
                picard_path = picard_path
            }

            call AlignAndSortBAM {
                input:
                sample_name = sample_name,
                fq1 = SamToFastq.fq1,
                fq2 = SamToFastq.fq2,

                ref = ref,
                sa = GenerateRefFiles.ref_sa,
                bwt = GenerateRefFiles.ref_bwt,
                amb = GenerateRefFiles.ref_amb,
                ann = GenerateRefFiles.ref_ann,
                pac = GenerateRefFiles.ref_pac,
                dict = GenerateRefFiles.ref_dict,
                fai = GenerateRefFiles.ref_index,
                docker = docker,
                mem_size_gb = small_mem_size_gb,
                disk_size = large_disk_size,
                picard_path = picard_path
            }
        }

        call MarkDuplicates {
            input:
            sample_name = sample_name,
            #sorted_bam = select_first([AlignAndSortBAM.bam, input_bam]),
            sorted_bam = input_bam,
            docker = docker,
            picard_path = picard_path,
            mem_size_gb = small_mem_size_gb,
            disk_size = disk_size
        }

        call ReorderBam {
            input:
            bam = MarkDuplicates.bam,
            ref = GenerateRefFiles.reference_fasta,
            dict = GenerateRefFiles.ref_dict,
            docker = docker,
            picard_path = picard_path,
            mem_size_gb = med_mem_size_gb,
            disk_size = disk_size
        }

        call HaplotypeCaller {
            input:
            input_bam = ReorderBam.out,
            input_bam_index = ReorderBam.out_index,
            sample_name = sample_name,
            gvcf_name = "${sample_name}.g.vcf",
            gvcf_index = "${sample_name}.g.vcf.idx",
            ref = GenerateRefFiles.reference_fasta,
            ref_dict = GenerateRefFiles.ref_dict,
            ref_index = GenerateRefFiles.ref_index,
            mem_size_gb = med_mem_size_gb,
            disk_size = med_disk_size,
            docker = docker,
            gatk_path = gatk_path
        }
    }

    call CombineGVCFs {
        input:
        vcf_files = HaplotypeCaller.output_gvcf,
        vcf_index_files = HaplotypeCaller.output_gvcf_index,
        ref = GenerateRefFiles.reference_fasta,
        ref_dict = GenerateRefFiles.ref_dict,
        ref_index = GenerateRefFiles.ref_index,
        docker = docker,
        gatk_path = gatk_path,
        mem_size_gb = med_mem_size_gb,
        disk_size = disk_size
    }

    call GenotypeGVCFs {
        input:
        vcf_file = CombineGVCFs.out,
        vcf_index_file = CombineGVCFs.out_index,
        ref = GenerateRefFiles.reference_fasta,
        ref_dict = GenerateRefFiles.ref_dict,
        ref_index = GenerateRefFiles.ref_index,
        docker = docker,
        gatk_path = gatk_path,
        mem_size_gb = med_mem_size_gb,
        disk_size = disk_size
    }

    call HardFiltration {
        input:
        vcf = GenotypeGVCFs.output_vcf_name,
        vcf_index = GenotypeGVCFs.output_vcf_index_name,
        snp_filter_expr = snp_filter_expr,
        indel_filter_expr = indel_filter_expr,
        ref = GenerateRefFiles.reference_fasta,
        ref_dict = GenerateRefFiles.ref_dict,
        ref_index = GenerateRefFiles.ref_index,
        output_filename = "${run_name}.hard_filtered.vcf.gz",
        docker = docker,
        gatk_path = gatk_path,
        mem_size_gb = extra_large_mem_size_gb,
        disk_size = extra_large_disk_size
    }

    call VcfToMSA {
        input:
            vcf = HardFiltration.out,
            ref = GenerateRefFiles.reference_fasta,
            output_filename = "${run_name}.msa.fasta",
            disk_size = disk_size,
            mem_size_gb = med_mem_size_gb
    }
    call IqTree2 {
        input:
        alignment = VcfToMSA.alignment,
        cluster_name = run_name,
        iqtree2_model = iqtree2_model,
        iqtree2_bootstraps = iqtree2_bootstraps,
        alrt = alrt,
        iqtree2_opts = iqtree2_opts,
        disk_size = disk_size,
        cpu = 4,
        memory = med_mem_size_gb
    }

    output {
        File gvcf = HardFiltration.out
    }
}

## TASK DEFINITIONS
task GbffToFasta {
    File ref_gbff
    String ref_gbff_basename = basename(ref_gbff, ".fasta")

    Int disk_size = 50
    Int mem_size_gb = 16
    String docker = "python:3.11-slim"

    command <<<
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
task SamToFastq {
    File in_bam
    String sample_name

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    String out_fq1 = "${sample_name}.1.fq"
    String out_fq2 = "${sample_name}.2.fq"
    command {
        java -Xmx${mem_size_gb}G -jar ${picard_path} SamToFastq INPUT=${in_bam} FASTQ=${out_fq1} SECOND_END_FASTQ=${out_fq2} VALIDATION_STRINGENCY=LENIENT
    }
    output {
        String done = "Done"
        File fq1 = out_fq1
        File fq2 = out_fq2
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        in_bam: "The bam file to convert to fastq."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        out_fq1: "The fastq file containing the first read of each pair."
        out_fq2: "The fastq file containing the second read of each pair"
    }
}


task AlignAndSortBAM {
    String sample_name

    File fq1
    File fq2

    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    String read_group = "'@RG\\tID:FLOWCELL_${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:LIB_${sample_name}'"
    command {
        bwa mem -R ${read_group} ${ref} ${fq1} ${fq2} | samtools view -bS -> ${sample_name}.aligned.bam
        java -Xmx${mem_size_gb}G -jar ${picard_path} SortSam I=${sample_name}.aligned.bam O=${sample_name}.sorted.bam SO=coordinate
    }

    output {
        File bam = "${sample_name}.sorted.bam"
    }

    runtime {
        task_name: "AlignBAM"
        preemptible: 5
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk "+ disk_size + " HDD"
    }

    parameter_meta {
        ref: "fasta file of reference genome"
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        fq_array: "An array containing the paths to the first and second fastq files."
        read_group: "The read group string that will be included in the bam header."
    }
}


# mark duplicate reads in bam
task MarkDuplicates {
    File sorted_bam
    String sample_name

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        set -euo pipefail

        java -Xmx${mem_size_gb}G -jar ${picard_path} MarkDuplicates \
            I=${sorted_bam} \
            O=${sample_name}.marked_duplicates.bam \
            M=${sample_name}.marked_duplicates.metrics
    }

    output {
        File bam = "${sample_name}.marked_duplicates.bam"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# reorder and index a bam
task ReorderBam {
    File ref
    File dict
    File bam
    String bam_prefix = basename(bam, '.bam')

    Int disk_size
    Int mem_size_gb
    String docker
    String picard_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # reorder bam
        java -Xmx${cmd_mem_size_gb}G -jar ${picard_path} ReorderSam \
            I=${bam} \
            O=${bam_prefix}.reordered.bam \
            R=${ref}

        # then index
        java -Xmx${cmd_mem_size_gb}G -jar ${picard_path} BuildBamIndex \
            I=${bam_prefix}.reordered.bam
    }

    output {
        File out = "${bam_prefix}.reordered.bam"
        File out_index = "${bam_prefix}.reordered.bai"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# merge vcfs before genotyping
task CombineGVCFs {
    File ref
    File ref_dict
    File ref_index
    Array[File] vcf_files
    Array[File] vcf_index_files

    String gvcf_out = "combined_gvcfs.vcf.gz"
    String gvcf_out_index = "combined_gvcfs.vcf.gz.tbi"

    Int disk_size
    Int mem_size_gb
    String docker
    String gatk_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T CombineGVCFs \
            -R ${ref} \
            -o ${gvcf_out} \
            --variant ${sep=" --variant " vcf_files}
    }
    output {
       File out = gvcf_out
       File out_index = gvcf_out_index
    }

    runtime {
        preemptible: 4
        docker:docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# genotype gvcfs
task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
    File vcf_file
    File vcf_index_file
    String gvcf_out = "genotyped_gvcfs.vcf.gz"

    Int disk_size
    Int mem_size_gb
    String docker
    String gatk_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T GenotypeGVCFs \
            -R ${ref} \
            -o ${gvcf_out} \
            --variant ${vcf_file} \
    }
    output {
        File output_vcf_name = gvcf_out
        File output_vcf_index_name = "${gvcf_out}.tbi"
    }

    runtime {
        preemptible: 4
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


task HardFiltration {
    # hard-filter a vcf, if vqsr not available
    # http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
    File ref
    File vcf
    File vcf_index
    File ref_dict
    File ref_index
    String output_filename

    String snp_filter_expr
    String indel_filter_expr

    Int disk_size
    Int mem_size_gb
    String docker
    String gatk_path

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # select snps
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path}\
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.g.vcf

        # filter indels
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path}\
            -T CombineVariants \
            -R ${ref} \
            --variant filtered_snps.g.vcf \
            --variant filtered_indels.g.vcf \
            -o ${output_filename} \
            --genotypemergeoption UNSORTED
    }

    output {
        File out = "${output_filename}"
    }

    runtime {
        preemptible: 3
        docker:docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# run haplotype caller
task HaplotypeCaller {
  File input_bam
  File input_bam_index

  String gvcf_name
  String gvcf_index

  File ref_dict
  File ref
  File ref_index
  String sample_name

  Int disk_size
  Int mem_size_gb
  Int cmd_mem_size_gb = mem_size_gb - 1

  String docker
  String gatk_path

  String out = "${sample_name}.g.vcf"

  command {
    java -Xmx${cmd_mem_size_gb}G -jar ${gatk_path} \
      -T HaplotypeCaller \
      -R ${ref} \
      -I ${input_bam} \
      -o ${gvcf_name} \
      -ERC "GVCF" \
      -ploidy 1 \
      -variant_index_type LINEAR \
      -variant_index_parameter 128000 \
      --read_filter OverclippedRead
  }

  output {
      #To track additional outputs from your task, please manually add them below
      File output_gvcf = "${gvcf_name}"
      File output_gvcf_index = "${gvcf_index}"
  }

  runtime {
    task_name: "HaplotypeCaller"
    preemptible: 3
    docker: docker
    memory: mem_size_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  parameter_meta {
      gatk: "Executable jar for the GenomeAnalysisTK"
      ref: "fasta file of reference genome"
      sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
      sample_dir: "The sample-specific directory inside output_dir for each sample."
      in_bam: "The bam file to call HaplotypeCaller on."
      out: "VCF file produced by haplotype caller."
  }
}

task VcfToMSA {
    File vcf
    File ref
    String vcf_basename = basename(vcf, ".vcf.gz")
    String output_filename

    Int disk_size
    Int mem_size_gb

    String docker = "biocontainers/bcftools:v1.9-1-deb_cv1"

    command <<<

        set -euo pipefail

        # index the VCF
        cp ${vcf} ${vcf_basename}.vcf.gz
        ls -lh
        echo "indexing the vcf"
        bcftools index -t ${vcf_basename}.vcf.gz
        ls -lh

        # Convert VCF to FASTA alignment
        bcftools query -l ${vcf} > sample_list.txt

        # Create consensus sequence for each sample
        while read sample; do
            # Create a consensus sequence for this sample
            bcftools consensus -f ${ref} -s $sample ${vcf} > $sample.fasta

            # Format header to just sample name
            sed -i "s/>.*/>$sample/" $sample.fasta
        done < sample_list.txt

        # Combine all individual FASTA files
        cat *.fasta > ${output_filename}
    >>>

    output {
        File alignment = "${output_filename}"
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
    Int disk_size = 100
    Int cpu = 4
    Int memory = 32

  command <<<
    # date and version control
    date | tee DATE

    # multiple sed statements to get down to a string that is just "version 2.1.2"
    iqtree2 --version | grep version | sed 's|.*version|version|;s| COVID-edition for Linux.*||' | tee VERSION

    # check if iqtree2_model input is set and output for sanity
    if [ -n "~{iqtree2_model}" ]; then
      echo "DEBUG: User provided iqtree2_model ~{iqtree2_model}, will use this for running iqtree2"
      IQTREE2_MODEL="~{iqtree2_model}"
    else
      echo "DEBUG: User did not supply an iqtree2_model input, will use iqtree2's model finder"
    fi

    # sanity check
    echo "DEBUG: IQTREE2_MODEL is set to: " $IQTREE2_MODEL

    # make sure there are more than 3 genomes in the dataset
    numGenomes=$(grep -o '>' ~{alignment} | wc -l)
    if [ "$numGenomes" -gt 3 ]; then
      cp ~{alignment} ./msa.fasta

      # run iqtree2
      #   -nt : number of CPU cores for multicore version
      #   -s : input alignment file
      #   -m : model
      #   -bb : number of bootstrap replicates
      #   -alrt : number of replicates to perform SH-like approximate likelihood ration test
      if [[ -v IQTREE2_MODEL ]] ; then # iqtree2 model set; use -m tag
        echo "DEBUG: running iqtree2 with the -m flag which is used to provide a model; user-specified " $IQTREE2_MODEL
        iqtree2 \
          -nt AUTO \
          -s msa.fasta \
          -m $IQTREE2_MODEL \
          -bb ~{iqtree2_bootstraps} \
          -alrt ~{alrt} ~{iqtree2_opts}

        # write the iqtree2_model used to a txt file for output as a string
        echo $IQTREE2_MODEL | tee IQTREE2_MODEL.TXT

      else # iqtree model is not set; do not use -m tag
        echo "DEBUG: running iqtree2 without the -m flag which is used to provide a model. Will default to iqtree2 default (Model Finder)"
        iqtree2 \
          -nt AUTO \
          -s msa.fasta \
          -bb ~{iqtree2_bootstraps} \
          -alrt ~{alrt} ~{iqtree2_opts}

        # for scenario where user did not specify iqtree2_model input nor core_genome boolean input, determine iqtree2_model used by parsing log file
        # first sed is to remove "Best-fit model: " and second sed is to remove anything after the word "chosen *", leaving only the name of the model
        grep "Best-fit model" msa.fasta.log | sed 's|Best-fit model: ||g;s|chosen.*||' | tee IQTREE2_MODEL.TXT

      fi

      # rename the final output newick file
      cp -v msa.fasta.contree ~{cluster_name}_iqtree.nwk
    else
      echo "DEBUG: not enough genomes provided; more than 3 are required to run iqtree2"
    fi
  >>>
  output {
    String date = read_string("DATE")
    String iqtree2_version = read_string("VERSION")
    File ml_tree = "~{cluster_name}_iqtree.nwk"
    String iqtree2_model_used = read_string("IQTREE2_MODEL.TXT")
    String iqtree2_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}