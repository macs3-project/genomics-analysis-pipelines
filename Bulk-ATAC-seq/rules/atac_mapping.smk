if mapper == "minimap2":
    rule atac_map:
        input:
            fastq1 = "%s/{fastqid}_R1.fastq.gz" % (config["fastqdir"]),
            fastq2 = "%s/{fastqid}_R2.fastq.gz" % (config["fastqdir"]),
        output:
            bam = temp("{OUT_DIR}/Alignment/{fastqid}.sortedByPos.bam")
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
            fastq1 = "%s/{fastqid}_R1.fastq.gz" % (config["fastqdir"]),
            fastq2 = "%s/{fastqid}_R2.fastq.gz" % (config["fastqdir"])
        output:
            bam = temp("{OUT_DIR}/Alignment/{fastqid}.sortedByPos.bam")
        params:
            genome = config["genome"]["bwaindex"],
        threads:
            config["options"]["cores"]
        shell:
            "bwa mem -t {threads} {params.genome} {input.fastq1} {input.fastq2} "
            "| samtools view --threads {threads} -b"
            "| samtools sort --threads {threads} -o {output.bam}"

