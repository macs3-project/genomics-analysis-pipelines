# mark duplicated reads with picard
rule atac_bammkdp:
    input:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.bam"
    output:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.mkdp.bam",
        metric = "{OUT_DIR}/Alignment/{name}.sortedByPos.mkdp.txt",
        tmp = temp(directory("{OUT_DIR}/Tmp/{name}"))
    shell:
        "picard MarkDuplicates -INPUT {input.bam} -OUTPUT {output.bam} -METRICS_FILE {output.metric} -TMP_DIR {output.tmp};"
        "rm {input.bam}"

# generate clean (filtered/deduplicated and chrM-removed and Q30 filtered ) bam file for macs
rule atac_bamrmdp:
    input:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.mkdp.bam",
    output:
        bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.bam",
    params:
        chrombed = config["annotation"]["chromBed"],
    shell:
        "samtools view --threads {threads} -b -L {params.chrombed} -F 3340 -f 0x2 -q 30 -o {output.bam} {input.bam};"

