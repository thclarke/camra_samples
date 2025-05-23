version 1.0

workflow combine_files_from_rows {
  input {
    Array[File] input_files  # Files from multiple table rows, passed as inputs
    String output_prefix     # Optional name to prefix the output
  }

  call collect_files {
    input: files = input_files,
           prefix = output_prefix
  }

  output {
    Array[File] combined_files = collect_files.collected_files
  }
}

task collect_files {
  input {
    Array[File] files
    String prefix
  }

  command {
    mkdir -p output
    for f in ~{sep=' ' files}; do
      cp "$f" output/
    done
  }

  output {
    Array[File] collected_files = glob("output/*")
  }

  runtime {
    docker: "ubuntu"
  }
}