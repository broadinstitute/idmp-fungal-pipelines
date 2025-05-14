version 1.0

task VCFToFasta {
    input {
    File vcf_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")

    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/vcftomsa:1.0.0"
    Int cpu = 4

    Int disk_size_gb = ceil((size(vcf_file, "GiB")) * 2) + 20
    Int mem_size_gb = ceil(size(vcf_file, "GiB") * 2.5) + 10

    }
    command <<<
    python3 /app/vcf2matrix.py \
      -f \
      -i ~{vcf_file} \
      --output-prefix ~{vcf_basename}
    >>>
    output {
        File alignment_fasta = "~{vcf_basename}.min4.fasta"
    }
    runtime {
        docker: docker
        memory: mem_size_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        disk: disk_size_gb + " GB"
        preemptible: 0
    }
}
