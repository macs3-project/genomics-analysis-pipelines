# wild card is "i"

if mapper == "minimap2":
    rule atac_map:
        input:
            fastq1 = "%s/%s_R1.fastq.gz" % (config["fastqdir"], SAMPLE1_FILES[{i}-1]),
            fastq2 = "%s/%s_R2.fastq.gz" % (config["fastqdir"], SAMPLE1_FILES[{i}-1]),
        output:
            bam = temp("{OUT_DIR}/Alignment/%s.sortedByPos.bam" % SAMPLE1_NAMES[{i}-1])
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
            fastq1 = "%s/%s_R1.fastq.gz" % (config["fastqdir"], SAMPLE1_FILES[{i}-1]),
            fastq2 = "%s/%s_R2.fastq.gz" % (config["fastqdir"], SAMPLE1_FILES[{i}-1]),
        output:
            bam = temp("{OUT_DIR}/Alignment/%s.sortedByPos.bam" % SAMPLE1_NAMES[{i}-1])
        params:
            genome = config["genome"]["bwaindex"],
        threads:
            config["options"]["cores"]
        shell:
            "bwa mem -t {threads} {params.genome} {input.fastq1} {input.fastq2} "
            "| samtools view --threads {threads} -b"
            "| samtools sort --threads {threads} -o {output.bam}"

