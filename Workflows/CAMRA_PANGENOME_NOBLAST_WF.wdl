version 1.0


task run_Pangenome {
    input {
       Array[File] gb_files
       File combined_blast
       Int hdd_sz
       Int gb_req
    }
    runtime{
        docker: "thclarke/pangenomepipeline:latest"
        memory: "~{gb_req} GB"   # Request 4 GB of memory
        disks: "local-disk ~{hdd_sz} SSD"

    }
    command <<< 
        mkdir ./gb_dir
        touch gb.list
        echo ~{sep = " " gb_files}
        for fl in ~{sep = " " gb_files}; do
          echo $fl
          cp $fl ./gb_dir/
          echo $(pwd)
          echo $(basename $fl) | awk -F "\." -v dir="$(pwd)" '{ print($1"\t"dir"/"$0); }' >> gb.list
          echo $(basename $fl) | awk -F "\." -v dir="$(pwd)" '{ echo $1"\t"dir"/"$0; }'
        done
        echo ~{combined_blast} 
        cp ~{combined_blast} ./
        perl /pangenome/bin/run_pangenome.pl --gb_dir ./gb_dir/ --no_blast --no_grid --panoct_local
    >>>
    output {
        File pangenome_centroids = "results/centroids.fasta.cleaned"
        File pangenome_stats = "results/overview_stats.txt"
        File pangenome_presence = "results/all_clusters_member_presence.txt"
        File pangenome_fs ="results/frameshifts.txt"
        File pangenome_roles = "results/centroids.cluster_roles.txt"
        File pangenome_aro = "results/aro_searches/centroids.fasta.cleaned.aro.txt"
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

    input {
        Array[File] gb_files
        String db_type
        File blast_output
        Int hdd_sz
	      Int gb_req
    }

    call run_Pangenome {
        input:
            gb_files = gb_files,
            combined_blast = blast_output,
            gb_req = gb_req,
            hdd_sz = hdd_sz
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
