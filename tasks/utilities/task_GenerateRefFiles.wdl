version 1.0

task GenerateRefFiles {
    input {
    File ref_fasta
    File input_bam
    String ref_fasta_basename = basename(ref_fasta, ".fasta")

    Int disk_size = 50
    Int mem_size_gb = 16
    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/samtools-picard-bwa:1.0.0-0.7.15-2.26.3-1634165082"
    }
    command <<<
       set -euo pipefail
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