version 1.0

import "../Tasks/task_plasmidfinder.wdl" as plasmidfinder
import "../Tasks/task_bvbrc.wdl" as bvbrc
import "../Tasks/task_pgap.wdl" as pgap

workflow annotation_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Annotation on an Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        String scientific_name # Example "Escherichia coli"

        File pgap_supplemental_data
        File local_assembly_path

        # Use if the assembly is stored in BV-BRC
      }   

    call pgap.pgap_annotate {
    input:
        sample_id = sample_name,
        assembly_fasta = local_assembly_path,
        supplemental_data = pgap_supplemental_data, # Refer to https://github.com/ncbi/pgap/wiki/Installation
        taxon = scientific_name, # Example "Escherichia coli"
        pgap_version = "2026-04-27.build8516"
    }
    output {

        File? pgap_ani_tax_report = pgap_annotate.pgap_ani_tax_report
        File? pgap_ani_tax_report_xml = pgap_annotate.pgap_ani_tax_report_xml
        File? pgap_annotated_proteins_fasta = pgap_annotate.pgap_annotated_proteins_fasta
        File? pgap_genomic_fasta = pgap_annotate.pgap_genomic_fasta
        File? pgap_gbk_annotation = pgap_annotate.pgap_gbk_annotation
        File? pgap_gff3_annotation = pgap_annotate.pgap_gff3_annotation
        File? pgap_annotated_sequence = pgap_annotate.pgap_annotated_sequence
        File? pgap_annotated_cds = pgap_annotate.pgap_annotated_cds
        File? pgap_translated_cds = pgap_annotate.pgap_translated_cds
        File? pgap_gff_with_sequence = pgap_annotate.pgap_gff_with_sequence
        File? pgap_foreign_sequence = pgap_annotate.pgap_foreign_sequence
        File? pgap_completeness_contamination = pgap_annotate.pgap_completeness_contamination
        File? pgap_log = pgap_annotate.pgap_log

    }
}
