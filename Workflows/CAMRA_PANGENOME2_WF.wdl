version 1.0

#import "../Tasks/task_pangenome.wdl" as pangenome



task makeFastaFiles {
  input {
    Array[File] gb_files
    String db_type
  }
  command <<<
     mkdir /tmp/gb_dir/
     mkdir /tmp/fasta_dir/

     for fl in ~{sep = " " gb_files}; do
        cp $fl /tmp/gb_dir/
        echo "/tmp/gb_dir/"$(basename $fl) >> /tmp/genomes.list
        echo "/tmp/gb_dir/"$(basename $fl)
     done
     cd /tmp/
     /pangenome/bin/parse_genbank_files.pl -l /tmp/genomes.list  -o  ./ --no_dos2unix
     for f in ./*pep; do
        cat $f >> /tmp/all_sequences.fasta
     done
     makeblastdb -in /tmp/all_sequences.fasta -dbtype ~{db_type} -out blast_db
     find /tmp/ -type f -name "blast_db*" > blast_list.txt
     find /tmp/fasta_dir/ -type f -name "GCF*pep" > fasta_list.txt



  >>>

  output {
	Array[File] blast_files = read_lines("blast_list.txt")
	String blast_db_prefix = "blast_db"
	Array[File] input_fastas = read_lines("fasta_list.txt")

  }

  runtime{
        docker: "thclarke/pangenomepipeline:latest"
	mem: "20 GiB"
  }
}



task RunBlast {
  input {
    File query_fasta
    String blast_db_prefix
    String db_type
    Array[File] blast_files
  }

  command <<<
    for fl in ~{sep = " " blast_files}; do
        cp $fl ./
        echo $fl
    done
    blastp -query ~{query_fasta} -db ~{blast_db_prefix} -out combined.blast -outfmt 6
  >>>

  output {
    File blast_output = "combined.blast"
  }

  runtime {
    docker: "thclarke/pangenomepipeline:latest"
  }
}

task makeCombinedBlast {
  input {
    Array[File]  blast_output
  }

  command <<<
     cat ~{ sep = ' '  blast_output} > combined.blast
  >>>

  output {
    File blast_output = "combined.blast"
  }

  runtime {
    docker: "thclarke/pangenomepipeline:latest"
  }
}


task run_Pangenome {
    input {
       Array[File] gb_files
       File combined_blast
    }
    runtime{
        docker: "thclarke/pangenomepipeline:latest"
         #disks: "local-disk 100 SSD"

    }
   	command <<<
        mkdir /tmp/gb_dir/
        for fl in ~{sep = " " gb_files}; do
          cp $fl /tmp/gb_dir/
        done
        cp ~{combined_blast} ./
        /pangenome/bin/run_pangenome.pl --gb_dir /tmp/gb_dir/ --no_blast --no_grid --panoct_local
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
    }


    call makeFastaFiles{
     input:
        gb_files = gb_files,
        db_type = db_type
  }


  scatter(fasta_file in makeFastaFiles.input_fastas) {
    call RunBlast {
      input:
        query_fasta = fasta_file,
        blast_db_prefix = makeFastaFiles.blast_db_prefix,
        blast_files = makeFastaFiles.blast_files,

        db_type = db_type
      }
    }

    call makeCombinedBlast {
        input:
          blast_output = RunBlast.blast_output

    }

    call run_Pangenome {
        input:
            gb_files = gb_files,
            combined_blast = makeCombinedBlast.blast_output
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
