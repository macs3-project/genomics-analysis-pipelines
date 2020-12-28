rule atac_qcstat:
    input:
        dirty_bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.mkdp.bam",
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        qc_stat = "{OUT_DIR}/QC/{name}.stat.txt",
        uniq_bam = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.unique.bam",
        uniq_bed = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.unique.bed",
        uniq_clean_bed = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.unique.clean.bed",
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
        insertl = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.unique.bam.insertl",
        insertlsum = "{OUT_DIR}/Alignment/{name}.sortedByPos.rmdp.clean.unique.bam.insertl.txt",
    shell:
        "samtools view {input.clean_bam} | cut -f 9 | awk '$1>0{{print}}' > {output.insertl};"
        "perl -e 'while(<>){{chomp;$bin=int($_/10);$count{{$bin}}+=1}}foreach my $key (sort {{$a<=>$b}} keys %count){{print $key*10,\"\\t\",$count{{$key}},\"\\n\"}}' {output.insertl} > {output.insertlsum};"
