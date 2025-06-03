version 1.0

task resolve_samn_to_assembly {
  input {
    String samn_id
  }

  command <<<
    esearch -db biosample -query ~{samn_id} |
    elink -target assembly |
     efetch -format docsum  > out.txt
     cat out.txt | xtract -pattern DocumentSummary -element AssemblyAccession > assembly.txt

      cat out.txt |xtract -pattern DocumentSummary -element FtpPath_RefSeq > genbank.txt
     if (du genbank.txt | cut -f1 == '0'); do 
         cat out.txt |xtract -pattern DocumentSummary -element FtpPath_GenBank > genbank.txt
     done;
     cat out.txt | xtract -pattern DocumentSummary -element AssemblyName > AssemblyName.txt
    >>>



  output {
    String assembly_id = read_string("assembly.txt")
    String genbank_path = read_string("genbank.txt")
    String AssemblyName = read_string("AssemblyName.txt")

  }

  runtime {
    docker: "ncbi/edirect"
  }
}

task download_genbank {
  input {
    String assembly_id
    String genbank_path
    String AssemblyName
  }

  command <<<
    python3 <<CODE
    import os
    import urllib.request

    prefix = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    assembly_id = "~{assembly_id}"
    genbank_path = "~{genbank_path}"
    AssemblyName = "~{AssemblyName}"
    parts = assembly_id.split("_")[1]
    dir_structure = "/".join([parts[i:i+3] for i in range(0, len(parts)-1, 3)])
    url = f"{genbank_path}/{assembly_id}_{AssemblyName}_genomic.gbff.gz"
    print(url)
    urllib.request.urlretrieve(url, f"{assembly_id}.gbff.gz")
    CODE
    gunzip ~{assembly_id}.gbff.gz
  >>>

  output {
    File genbank_file = "~{assembly_id}.gbff"
  }

  runtime {
    docker: "python:3"
  }
}

workflow fetch_genbank {

meta {
        author: "Thomas Clarke"
        email: "tclarke@jcvi.org"
        description: "Run Pangenome."
  }

  input {
    String samn_id
  }

  call resolve_samn_to_assembly {
    input:
      samn_id = samn_id
  }

  call download_genbank {
    input:
      assembly_id = resolve_samn_to_assembly.assembly_id,
      genbank_path = resolve_samn_to_assembly.genbank_path,
      AssemblyName = resolve_samn_to_assembly.AssemblyName
  }

  output {
   File  genbank_file = download_genbank.genbank_file
  }

}


