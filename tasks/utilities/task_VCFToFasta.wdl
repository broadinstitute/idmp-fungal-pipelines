task VCFToFasta {
    File vcf_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")

    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/vcftomsa:1.0.0"
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
        docker: docker
        memory: memory + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB"
        preemptible: 0
    }
}
