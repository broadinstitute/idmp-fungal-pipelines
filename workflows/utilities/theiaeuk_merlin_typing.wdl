version 1.0

# theiaeuk
import "../../tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl" as snippy_gene_query
import "../../tasks/gene_typing/variant_detection/task_snippy_variants.wdl" as snippy
import "../../tasks/species_typing/candida/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/quality_control/read_filtering/task_filter_bam.wdl" as bam_filter

workflow theiaeuk_merlin_typing {
    meta {
        description: "Workflow for bacterial and fungal species typing; based on the Bactopia subworkflow Merlin (https://bactopia.github.io/bactopia-tools/merlin/)"
    }
    input {
        String samplename
        String merlin_tag
        File assembly
        File? read1
        File? read2
        # subworkflow logic
        Boolean assembly_only = false
        Boolean ont_data = false
        Boolean theiaeuk = false
        # docker options
        String? cauris_cladetyper_docker_image
        String? snippy_gene_query_docker_image
        String? snippy_variants_docker_image
        # cladetyper options - primarily files we host
        Int? cladetyper_kmer_size
        File? cladetyper_ref_clade1
        File? cladetyper_ref_clade1_annotated
        File? cladetyper_ref_clade2
        File? cladetyper_ref_clade2_annotated
        File? cladetyper_ref_clade3
        File? cladetyper_ref_clade3_annotated
        File? cladetyper_ref_clade4
        File? cladetyper_ref_clade4_annotated
        File? cladetyper_ref_clade5
        File? cladetyper_ref_clade5_annotated
        # snippy options - mostly files we host
        String? snippy_query_gene
        Int? snippy_map_qual
        Int? snippy_base_quality
        Int? snippy_min_coverage
        Float? snippy_min_frac
        Int? snippy_min_quality
        Int? snippy_maxsoft
    }
    # theiaeuk
    if (theiaeuk) {
        if (merlin_tag == "Candidozyma auris" || merlin_tag == "Candida auris") {
            call cauris_cladetyper.cauris_cladetyper as cladetyper {
                input:
                    assembly_fasta = assembly,
                    samplename = samplename,
                    kmer_size = cladetyper_kmer_size,
                    ref_clade1 = cladetyper_ref_clade1,
                    ref_clade1_annotated = cladetyper_ref_clade1_annotated,
                    ref_clade2 = cladetyper_ref_clade2,
                    ref_clade2_annotated = cladetyper_ref_clade2_annotated,
                    ref_clade3 = cladetyper_ref_clade3,
                    ref_clade3_annotated = cladetyper_ref_clade3_annotated,
                    ref_clade4 = cladetyper_ref_clade4,
                    ref_clade4_annotated = cladetyper_ref_clade4_annotated,
                    ref_clade5 = cladetyper_ref_clade5,
                    ref_clade5_annotated = cladetyper_ref_clade5_annotated,
                    docker = cauris_cladetyper_docker_image
            }
            if (!assembly_only && !ont_data) {
                call snippy.snippy_variants as snippy_cauris { # no ONT support right now
                    input:
                        reference_genome_file = cladetyper.annotated_reference,
                        read1 = select_first([read1]),
                        read2 = read2,
                        samplename = samplename,
                        map_qual = snippy_map_qual,
                        base_quality = snippy_base_quality,
                        min_coverage = snippy_min_coverage,
                        min_frac = snippy_min_frac,
                        min_quality = snippy_min_quality,
                        maxsoft = snippy_maxsoft,
                        docker = snippy_variants_docker_image
                }
                call snippy_gene_query.snippy_gene_query as snippy_gene_query_cauris {
                    input:
                        samplename = samplename,
                        snippy_variants_results = snippy_cauris.snippy_variants_results,
                        reference = cladetyper.annotated_reference,
                        query_gene = select_first([snippy_query_gene, "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,B9J08_005340,B9J08_000401,B9J08_003102,B9J08_003737,B9J08_005343"]),
                        docker = snippy_gene_query_docker_image
                }
                call bam_filter.bam_filter_fixmates as filter_bam {
                    input:
                        input_bam = snippy_cauris.snippy_variants_bam,
                        output_prefix = samplename
                }
            }
        }
    }
    output {
        # theiaeuk
        # c auris
        String? clade_type = cladetyper.gambit_cladetype
        String? cladetyper_version = cladetyper.gambit_version
        String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
        String? cladetype_annotated_ref = cladetyper.annotated_reference
        File? claderef_fasta = cladetyper.claderef_fasta
        # snippy variants
        String snippy_variants_query = select_first([snippy_gene_query_cauris.snippy_variants_query, "No matching taxon detected"])
        String snippy_variants_query_check = select_first([snippy_gene_query_cauris.snippy_variants_query_check, "No matching taxon detected"])
        String snippy_variants_hits = select_first([snippy_gene_query_cauris.snippy_variants_hits, "No matching taxon detected"])
        String snippy_variants_gene_query_results = select_first([snippy_gene_query_cauris.snippy_variants_gene_query_results, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
        File? filtered_bam = filter_bam.filtered_bam
    }
}