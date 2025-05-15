version 1.0


task run_Pangenome {
    input {
        String gb_dir
    }
    runtime{
        docker:'thclarke/pangenomepipeline:latest'
        disks:'/Users/thclarke/Downloads/PanGenome_Cromwell/'
    }
    command <<<
        mkdir /tmp/gb_direct/
        ls ~{gb_dir} > file_list.txt
        ls ~{gb_dir} | awk '{ system("cp ~{gb_dir}/"$0" /tmp/gb_direct/"); }'
        ls ~{gb_dir} | awk '{ system("echo ~{gb_dir}/"$0); }'

        #cp -r ${glob(gb_dir + "/*")} "/tmp/gb_direct/"
        /pangenome/bin/run_pangenome.pl --gb_dir /tmp/gb_direct/ --working_dir /tmp/ --blast_local --no_grid --panoct_local --working_dir .
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