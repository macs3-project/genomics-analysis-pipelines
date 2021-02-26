rule atac_qcstat:
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
    threads:
        config["options"]["cores"]
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_BulkQCStat.benchmark"
    shell:
        "echo 'flagstat:' > {output.qc_stat};"
        "samtools flagstat --threads {threads} {input.dirty_bam} >> {output.qc_stat};"
        "samtools view --threads {threads} -b -F 3340 -f 0x2 -q 30 -o {output.uniq_bam} {input.dirty_bam};"
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

rule atac_peakqc:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
        peakxls = "{OUT_DIR}/Analysis/{name}_peaks.xls",
    output:
        peak_qc = "{OUT_DIR}/QC/{name}.peakstat.txt",
    params:
        promoter = config["annotation"]["promoter"],
        chrMregion = config["annotation"]["MtBed"],
        blacklist = config["annotation"]["blacklist"],
        DHS = config["annotation"]["DHS"],
    threads:
        config["options"]["cores"],
    benchmark:
        "{OUT_DIR}/Benchmark/{name}_PeakQCStat.benchmark",
    shell:
        "grep 'total fragments in treatment' {input.peakxls} | perl -pe 's/# //' > {output.peak_qc};"
        "echo 'total number of peaks:' >> {output.peak_qc};"
        "wc -l {input.peak} | cut -f 1 -d' ' >> {output.peak_qc};"
        "echo 'number of peaks over FC 2:' >> {output.peak_qc};"
        "awk '$7>=2{{print}}' {input.peak} | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in blacklist regions:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.blacklist} -u | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in chrM:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.chrMregion} -u | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in promoter regions:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.promoter} -u | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in DHS regions:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.DHS} -u | wc -l >> {output.peak_qc};"

rule atac_frag:
    input:
        clean_bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.bam"
    output:
        insertl = temp("{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.unique.bam.insertl"),
        insertlsum = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.unique.bam.insertl.txt",
    shell:
        "samtools view {input.clean_bam} | cut -f 9 | awk '$1>0{{print}}' > {output.insertl};"
        "perl -e 'while(<>){{chomp;$bin=int($_/10);$count{{$bin}}+=1}}foreach my $key (sort {{$a<=>$b}} keys %count){{print $key*10,\"\\t\",$count{{$key}},\"\\n\"}}' {output.insertl} > {output.insertlsum};"

rule atac_plotfrag:
    input:
        fragstat1 = FRAG_STAT1,
        fragstat2 = FRAG_STAT2,
    output:
        png = FRAG_PNG,
    params:
        fragstat1param = ",".join(FRAG_STAT1),
        fragstat2param = ",".join(FRAG_STAT2),
        s1name = config["sample1"],
        s2name = config["sample2"],
    shell:
        "utils/fragstat.R -a {params.fragstat1param} -b {params.fragstat2param} -i {params.s1name} -j {params.s2name} -o {output.png};"

rule atac_sum_qcstat:
    input:
        qcstat1 = SEQ_STAT1,
        qcstat2 = SEQ_STAT2,
    output:
        seqqcsummary = SEQ_QC_SUMMARY,
    params:
        filelist = " ".join(SEQ_STAT1 + SEQ_STAT2)
    shell:
        "utils/atac_seqqc_summary.py {params.filelist} > {output.seqqcsummary};"

rule atac_sum_peakstat:
    input:
        peakstat1 = PEAK_STAT1,
        peakstat2 = PEAK_STAT2,
    output:
        peakqcsummary = PEAK_QC_SUMMARY,
    params:
        filelist = " ".join(PEAK_STAT1 + PEAK_STAT2),
    shell:
        "utils/atac_peakqc_summary.py {params.filelist} > {output.peakqcsummary};"

rule atac_plot_gss:
    input:
        bw1 = BIGWIG_SPMR1,
        bw2 = BIGWIG_SPMR2,
    output:
        mat = GSSMAT,
        heatmap = GSSHEATMAP,
        profile = GSSPROFILE,
    params:
        gtf = config["annotation"]["geneGTF"],
        bwlist = " ".join( BIGWIG_SPMR1 + BIGWIG_SPMR2 ),
    shell:
        "awk '$3==\"gene\"{{print}}' {params.gtf} | perl -ne 'chomp;@F=split(/\\t/);$g=$F[8];$g=~s/^gene_id\\ \\\"(\\S+)\\\".*/$1/;$c=$F[0];$s=$F[6];if ($s==\"+\"){{$p=$F[3]-1}}else{{$p=$F[4]}}print join(\"\\t\",$c,$p,$p+1,$g,\".\",$s),\"\\n\"' > gss.bed;"
        "computeMatrix reference-point -S {params.bwlist} -R gss.bed --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --skipZeros -o {output.mat};"
        "rm -f gss.bed;"
        "plotHeatmap -m {output.mat} -out {output.heatmap};"
        "plotProfile -m {output.mat} -out {output.profile} --plotType=se --perGroup;"
        
#rule atac_profile_gss:
#    input:
#        peak   = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
#        bigwig = "{OUT_DIR}/Analysis/{name}_spmr.bw",
#    output:
#        fig = "{OUT_DIR}/QC/{name}_profile_gss.png",
#    params:
#        gtf = config["annotation"]["geneGTF"],
#    shell:
#        "awk '$3==\"gene\"{{print}}' {params.gtf} | perl -ne 'chomp;@F=split(/\\t/);$s=$F[8];$s=~s/^gene_id\\ \\\"\(\\S+\)\\\".*/$1/;print $F[0],$s,\"\\n\"'"
