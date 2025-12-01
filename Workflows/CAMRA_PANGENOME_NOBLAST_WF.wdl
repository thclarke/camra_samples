version 1.0


task run_Pangenome {
    input {
       Array[File] gb_files
       File combined_blast
       String gb_name
       Int hdd_sz
       Int gb_req
    }
    runtime{
        docker: "thclarke/pangenomepipeline:latest"
        memory: "~{gb_req} GB"   # Request 4 GB of memory
        disks: "local-disk ~{ hdd_sz } HDD"  # Request 4 GB of memory
    }
    command <<< 
        #cd /tmp/
        mkdir gb_dir
        touch gb.list
        echo ~{sep = " " gb_files}
        for fl in ~{sep = " " gb_files}; do
          echo $fl
          cp "$fl" ./
          echo $(basename $fl) | awk -F'.' -v dir="$(pwd)" '{ system("cp "$0" "dir"/gb_dir/"$1".~{gb_name}"); }'
          #echo $(basename $fl) | awk -F'\.~{gb_name}' -v dir="$(pwd)" '{ system("mv "$0" "dir"/"$0); }'
        done
        echo ~{combined_blast} 
        cp ~{combined_blast} ./
        for f in *.~{gb_name}; do
            echo "$f"
            chmod +rw "$f"
            mv "$f" ./gb_dir/
            #if [ -e ./gb_dir/$(basename $f) ]; then
            echo "--"$(basename $f)
            echo $(basename $f) | awk -F'.~{gb_name}' -v dir="$(pwd)"  '{ print $1"\t"dir"/gb_dir/"$0; }' >> gb.list
            echo $(basename $f) | awk -F'.~{gb_name}' -v dir="$(pwd)" '{ system("echo "$1"\t"dir"/gb_dir/"$0); }'
            #fi
        done
        /pangenome/bin/core_hmm_checker_prep.pl -g gb.list  -o ./
        touch combined.fasta
        touch combined.att
        echo -n > gb.list
        cat  gb.list
        mkdir att_dir
        DIR="$PWD"
        echo $DIR
        cd pep
        for f in *patt; do
          chmod +rw $f
          echo "$f"
          if [  -s "$f" ]; then
             #cat "$f" >> ../combined.att
             #mv "$f" ../att_dir/
             echo $(basename $f) | awk -F'.patt' -v dir="$DIR" '{ print $1"\t"dir"/gb_dir/"$1".~{gb_name}"; }' >> ../gb.list
           fi
           if [ ! -s "$f" ]; then
              echo $(basename $f) | awk -F'.patt' -v dir="$DIR" '{ system("rm ../gb_dir/"$1".~{gb_name}"); }'
           fi
        done
        cd ../
        cat gb.list
        #perl /pangenome/bin/run_pangenome.pl --gb_ gb.list --no_blast --no_grid --panoct_local --combined_fasta ./combined.fasta --combined_att ./combined.att 
        perl /pangenome/bin/run_pangenome.pl --gb_list_file gb.list --no_blast --no_grid --panoct_local --panoct_verbose
    >>>
    output {
        File gb_list = "gb.list"
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
    	File pangenome_core_att = "results/fGIs/Core.attfGI";            
        File pangenome_shared = "results/shared_clusters.txt";             
        File pangenome_conses_in = "results/fGIs/consensus.txt";          
        File pangenome_fgi_insert = "results/fGIs/FGI_inserts.details";    
        File pangenome_fgi_att = "results/fGIs/fGI.att";          
        File pangenome_comb_att = "combined.att";                     
        File pangenome_sing_cluster = "results/singletons_clusters.txt";
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
        String gb_name
        File blast_output
        Int hdd_sz
	      Int gb_req
    }

    call run_Pangenome {
        input:
            gb_files = gb_files,
            combined_blast = blast_output,
            gb_req = gb_req,
            hdd_sz = hdd_sz,
            gb_name = gb_name
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
        File gb_list = run_Pangenome.gb_list
   
    }

}
