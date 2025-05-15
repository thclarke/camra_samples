version 1.0

#import "../Tasks/task_pangenome.wdl" as pangenome

task list_files {
  input {
    String gb_dir
  }
  command <<<
    ls ~{gb_dir} > file_list.txt
  >>>
  output {
    Array[File] files_in_dir = glob(gb_dir + "/*")
    File file_listing = "file_list.txt"
  }
}


task run_Pangenome {
    input {
        String gb_dir
    }
    #runtime{
         #docker: "thclarke/pangenomepipeline:latest"
         #disks: "local-disk 100 SSD"

    #}
    command <<<
        export PATH=/Users/thclarke/Downloads/ncbi-blast-2.2.28+/bin/:/Users/thclarke/Downloads/hmmer/src/:/Users/thclarke/Downloads/EMBOSS-6.6.0/emboss/:$PATH
        /Users/thclarke/Downloads/Git/PanGenomePipeline/pangenome/bin/run_pangenome.pl --gb_dir ~{gb_dir} --blast_local --no_grid --panoct_local --working_dir .
    >>>
    output {
        File file_listing = "file_list.txt"
        File pangenome_centroids = "results/centroids.fasta"
        File pangenome_stats = "results/overview_stats.txt"
        File pangenome_presence = "results/all_clusters_member_presence"
        File pangenome_fs ="results/frameshifts.txt"
        File pangenome_roles = "results/centroids.cluster_roles.txt"
        File pangenome_aro = "results/ centroids.fasta.cleaned.aro"
        File pangenome_fgi_index = "results/fGIs/fGI_report.txt.index"
        File pangenome_fgi_report = "results/fGIs/fGI_report.txt"
        File pangenome_tigrfam = "results/centroids.fasta_TIGRFAMs_15.0_HMM.tblout"
        File pangenome_pfam = "results/centroids.fasta_Pfam-A.tblout"
    }
}

workflow pangenome   {
    meta {
        author: "Thomas Clarke"
        email: "tclarke@jcvi.org"
        description: "Run Pangenome."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        String gb_dir
    }
    call list_files {
        input:
            gb_dir = gb_dir
    }


    call run_Pangenome {
        input:
            gb_dir = gb_dir
    }



    output {
        File pangenome_centorids = run_Pangenome.pangenome_centroids
        File pangenome_stats = run_Pangenome.pangenome_stats
        File pangenome_presence = run_Pangenome.pangenome_presence
        File pangenome_fs = run_Pangenome.pangenome_fs
        File pangenome_roles = run_Pangenome.pangenome_roles
        File pangenome_aro = run_Pangenome.pangenome_aro
        File pangenome_fgi_index = run_Pangenome.pangenome_fgi_index
        File pangenome_fgi_report = run_Pangenome.pangenome_fgi_report
        File pangenome_tigrfam = run_Pangenome.pangenome_tigrfam
        File pangenome_pfam = run_Pangenome.pangenome_pfam
    }

}