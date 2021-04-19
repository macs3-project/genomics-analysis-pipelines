# qc for sequencing alignment from each replicate
rule chip_qcstat:
    input:
        dirty_bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.mkdp.bam",
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        qc_stat = "{OUT_DIR}/QC/{name}.stat.txt",
        uniq_bam = temp("{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.unique.bam"),
        uniq_bed = temp("{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.unique.bed"),
        uniq_clean_bed = temp("{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.unique.clean.bed"),
    params:
        promoter = config["annotation"]["promoter"],
        chrMregion = config["annotation"]["MtBed"],
        bwflag_opt = bwflag_filter,
    threads:
        config["options"]["cores"],
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_BulkQCStat.benchmark"
    shell:
        "echo 'flagstat:' > {output.qc_stat};"
        "samtools flagstat --threads {threads} {input.dirty_bam} >> {output.qc_stat};"
        "samtools view --threads {threads} -b {params.bwflag_opt} -q 30 -o {output.uniq_bam} {input.dirty_bam};"
        "echo 'mapped Q30 reads:' >> {output.qc_stat};"
        "bedtools bamtobed -i {output.uniq_bam} > {output.uniq_bed};"
        "wc -l {output.uniq_bed} >> {output.qc_stat};"
        "echo 'chrM reads:' >> {output.qc_stat};"
        "bedtools intersect -a {output.uniq_bed} -b {params.chrMregion} -u | wc -l >> {output.qc_stat};"
        "echo 'non chrM reads:' >> {output.qc_stat};"
        "bedtools intersect -a {output.uniq_bed} -b {params.chrMregion} -v > {output.uniq_clean_bed};"
        "wc -l {output.uniq_clean_bed} >> {output.qc_stat};"
        "echo 'non chrM reads in promoter:' >> {output.qc_stat};"
        "bedtools intersect -a {output.uniq_clean_bed} -b {params.promoter} -u | wc -l >> {output.qc_stat};"
        "echo 'non chrM reads in peak:' >> {output.qc_stat};"
        "bedtools intersect -a {output.uniq_clean_bed} -b {input.peak} -u | wc -l >> {output.qc_stat};"

# qc for peakcalls from each replicate
def get_name ( wcs ):
    if config["options"]["paired"] == "Y":
        return "fragments"
    else:
        return "tags"
    
rule chip_peakqc:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
        peakxls = "{OUT_DIR}/Analysis/{name}_peaks.xls",
    output:
        peak_qc = "{OUT_DIR}/QC/{name}.peakstat.txt",
    params:
        name_tag =  get_name,
        promoter = config["annotation"]["promoter"],
        chrMregion = config["annotation"]["MtBed"],
        blacklist = config["annotation"]["blacklist"],
        DHS = config["annotation"]["DHS"],
    threads:
        config["options"]["cores"],
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_PeakQCStat.benchmark",
    shell:
        "grep 'total {params.name_tag} in treatment' {input.peakxls} | perl -pe 's/#\ //' > {output.peak_qc}; "
        "echo 'total number of peaks:' >> {output.peak_qc}; "
        "wc -l {input.peak} | cut -f 1 -d' ' >> {output.peak_qc}; "
        "echo 'number of peaks over FC 2:' >> {output.peak_qc}; "
        "awk '$7>=2{{print}}' {input.peak} | wc -l >> {output.peak_qc}; "
        "echo 'number of peaks in blacklist regions:' >> {output.peak_qc}; "
        "bedtools intersect -a {input.peak} -b {params.blacklist} -u | wc -l >> {output.peak_qc}; "
        "echo 'number of peaks in chrM:' >> {output.peak_qc}; "
        "bedtools intersect -a {input.peak} -b {params.chrMregion} -u | wc -l >> {output.peak_qc}; "
        "echo 'number of peaks in promoter regions:' >> {output.peak_qc}; "
        "bedtools intersect -a {input.peak} -b {params.promoter} -u | wc -l >> {output.peak_qc}; "
        "echo 'number of peaks in DHS regions:' >> {output.peak_qc}; "
        "bedtools intersect -a {input.peak} -b {params.DHS} -u | wc -l >> {output.peak_qc}; "

# combine all sequencing qc together
rule chip_sum_qcstat:
    input:
        qcstat_t = SEQ_STAT_T,
        qcstat_c = SEQ_STAT_C,
    output:
        seqqcsummary = SEQ_QC_SUMMARY,
    params:
        filelist = " ".join(SEQ_STAT_T + SEQ_STAT_C)
    shell:
        "utils/chip_seqqc_summary.py {params.filelist} > {output.seqqcsummary};"

# combine all peak qc together
rule chip_sum_peakstat:
    input:
        peakstat = PEAK_STAT,
    output:
        peakqcsummary = PEAK_QC_SUMMARY,
    params:
        filelist = " ".join(PEAK_STAT),
    shell:
        "utils/chip_peakqc_summary.py {params.filelist} > {output.peakqcsummary};"


rule chip_plot_gss:
    input:
        bw = BIGWIG_SPMR,
    output:
        mat = GSSMAT,
        heatmap = GSSHEATMAP,
        profile = GSSPROFILE,
    params:
        gtf = config["annotation"]["geneGTF"],
        bwlist = " ".join( BIGWIG_SPMR),
    shell:
        "awk '$3==\"gene\"{{print}}' {params.gtf} | perl -ne 'chomp;@F=split(/\\t/);$g=$F[8];$g=~s/^gene_id\\ \\\"(\\S+)\\\".*/$1/;$c=$F[0];$s=$F[6];if ($s==\"+\"){{$p=$F[3]-1}}else{{$p=$F[4]}}print join(\"\\t\",$c,$p,$p+1,$g,\".\",$s),\"\\n\"' > gss.bed;"
        "computeMatrix reference-point -S {params.bwlist} -R gss.bed --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros -o {output.mat};"
        "rm -f gss.bed;"
        "plotHeatmap -m {output.mat} -out {output.heatmap};"
        "plotProfile -m {output.mat} -out {output.profile} --plotType=se --perGroup;"
        
