if config["options"]["mapper"] == "minimap2":
    if config["options"]["paired"]:   
        rule chip_map:
            input:
                fastq1 = "{OUT_DIR}/Raw/{name}_R1.fastq.gz",
                fastq2 = "{OUT_DIR}/Raw/{name}_R2.fastq.gz",
            output:
                bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.bam"
            params:
                genome = config["genome"]["mmi"],
            threads:
                config["options"]["cores"]
            shell:
                "minimap2 -ax sr -t {threads} {params.genome} {input.fastq1} {input.fastq2} "
                "| samtools view --threads {threads} -b"
                "| samtools sort --threads {threads} -o {output.bam}"
    else:
        rule chip_map:
            input:
                fastq = "{OUT_DIR}/Raw/{name}.fastq.gz",
            output:
                bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.bam"
            params:
                genome = config["genome"]["mmi"],
            threads:
                config["options"]["cores"]
            shell:
                "minimap2 -ax sr -t {threads} {params.genome} {input.fastq}"
                "| samtools view --threads {threads} -b"
                "| samtools sort --threads {threads} -o {output.bam}"

elif config["options"]["mapper"] == "bwa-mem":
    if config["options"]["paired"]:
        rule chip_map:
            input:
                fastq1 = "{OUT_DIR}/Raw/{name}_R1.fastq.gz",
                fastq2 = "{OUT_DIR}/Raw/{name}_R2.fastq.gz",
            output:
                bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.bam"
            params:
                genome = config["genome"]["bwaindex"],
            threads:
                config["options"]["cores"]
            shell:
                "bwa mem -t {threads} {params.genome} {input.fastq1} {input.fastq2} "
                "| samtools view --threads {threads} -b"
                "| samtools sort --threads {threads} -o {output.bam}"
    else:
        rule chip_map:
            input:
                fastq = "{OUT_DIR}/Raw/{name}.fastq.gz",
            output:
                bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.bam"
            params:
                genome = config["genome"]["bwaindex"],
            threads:
                config["options"]["cores"]
            shell:
                "bwa mem -t {threads} {params.genome} {input.fastq} "
                "| samtools view --threads {threads} -b"
                "| samtools sort --threads {threads} -o {output.bam}"

