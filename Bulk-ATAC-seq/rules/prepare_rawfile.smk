rule prepare_link:
    input:
        fastq1 = "%s/%s_R1.fastq.gz" % (config["fastqdir"], sample_dict["{name}"]),
        fastq2 = "%s/%s_R2.fastq.gz" % (config["fastqdir"], sample_dict["{name}"]),
    output:
        fastqlink1 = "{OUT_DIR}/Raw/{name}_R1.fastq.gz",
        fastqlink2 = "{OUT_DIR}/Raw/{name}_R2.fastq.gz",
    shell:
        "ln -s {input.fastq1} {output.fastqlink1}; "
        "ln -s {input.fastq2} {output.fastqlink2}; "        
