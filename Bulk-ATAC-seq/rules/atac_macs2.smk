# mark duplicated reads with picard
rule atac_bammkdp:
    input:
        bam = "{OUT_DIR}/Alignment/{fastqid}.sortedByPos.bam"
    output:
        bam = "{OUT_DIR}/Alignment/{fastqid}.sortedByPos.mkdp.bam",
        metric = "{OUT_DIR}/Alignment/{fastqid}.sortedByPos.mkdp.txt",
        tmp = temp(directory("{OUT_DIR}/Tmp/{fastqid}"))
    shell:
        "picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR={output.tmp};"
        "rm {input.bam}"

# generate clean (filtered/deduplicated and chrM-removed and Q30 filtered ) bam file for macs
rule atac_bamrmdp:
    input:
        bam = "{OUT_DIR}/Alignment/{fastqid}.sortedByPos.mkdp.bam",
    output:
        bam = "{OUT_DIR}/Alignment/{fastqid}.sortedByPos.rmdp.clean.bam",
    params:
        chrombed = config["annotation"]["chromBed"],
    shell:
        "samtools view --threads {threads} -b -L {params.chrombed} -F 3340 -f 0x2 -q 30 -o {output.bam} {input.bam};"

# macs takes clean bam file to call peaks
rule atac_callpeak:
    input:
        bam = "{OUT_DIR}/Alignment/{fastqid}.sortedByPos.rmdp.clean.bam" 
    output:
        peak = "{OUT_DIR}/Analysis/{fastqid}_peaks.narrowPeak",
        bdg = "{OUT_DIR}/Analysis/{fastqid}_treat_pileup.bdg",
        xls = "{OUT_DIR}/Analysis/{fastqid}_peaks.xls",
    params:
        name = "{fastqid}",
    log:
        "{OUT_DIR}/Log/{fastqid}_macs2_peak.log"
    benchmark:
        "{OUT_DIR}/Benchmark/{fastqid}_callpeak.benchmark"
    shell:
        "macs2 callpeak -g hs --outdir {OUT_DIR}/Analysis -n {params.name} --keep-dup all -B -q 0.05 -f BAMPE --SPMR -t {input.bam};"

# bdg2bw converts bedGraph to bigWig files
rule atac_bdg2bw:
    input:
        bdg = "{OUT_DIR}/Analysis/{bdgname}.bdg",
    output:
        bw = "{OUT_DIR}/Analysis/{bdgname}.bw",
    params:
        chromlen = config["annotation"]["chromInfo"],
    shell:
        "bdg2bw {input.bdg} {params.chromlen};"