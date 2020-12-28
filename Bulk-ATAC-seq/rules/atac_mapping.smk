if mapper == "minimap2":
   rule atac_map:
        input:
            fastq1 = "{OUT_DIR}/Raw/{name}_R1.fastq.gz",
            fastq2 = "{OUT_DIR}/Raw/{name}_R2.fastq.gz",	
        output:
            bam = temp("{OUT_DIR}/Alignment/{name}.sortedByPos.bam")
        params:
            genome = config["genome"]["mmi"],
        threads:
            config["options"]["cores"]
        shell:
            "minimap2 -ax sr -t {threads} {params.genome} {input.fastq1} {input.fastq2} "
            "| samtools view --threads {threads} -b"
            "| samtools sort --threads {threads} -o {output.bam}"

elif mapper == "bwa-mem":
    rule atac_map:
        input:
            fastq1 = "{OUT_DIR}/Raw/{name}_R1.fastq.gz",
            fastq2 = "{OUT_DIR}/Raw/{name}_R2.fastq.gz",	
        output:
            bam = temp("{OUT_DIR}/Alignment/{name}.sortedByPos.bam")
        params:
            genome = config["genome"]["bwaindex"],
        threads:
            config["options"]["cores"]
        shell:
            "bwa mem -t {threads} {params.genome} {input.fastq1} {input.fastq2} "
            "| samtools view --threads {threads} -b"
            "| samtools sort --threads {threads} -o {output.bam}"

