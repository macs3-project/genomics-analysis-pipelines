# peak annotation with ChIPseeker
rule chip_peakanno:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        pdf = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak_chip_annotation.pdf"
    params:
        gmtfolder = config["annotation"]["GMT"],
    shell:
        "utils/chip_annotate.R {input.peak} {params.gmtfolder};"

# giggle search against CistromeDB
rule chip_giggle_cistrome:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        pdf = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak_giggle_cistrome.pdf"
    params:
        giggleexec = config["annotation"]["giggle"],
        gigglepath = config["annotation"]["giggledb"],
    shell:
        "utils/chip_giggle.R {input.peak} {params.giggleexec} {params.gigglepath};"

# homer
rule chip_homer:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        peaktop500 = temp("{OUT_DIR}/Analysis/{name}_peaks.narrowPeak.top500.bed"),
        homeroutput = temp("{OUT_DIR}/Analysis/{name}_HOMER"),
        homeroutputgz = "{OUT_DIR}/Analysis/{name}_HOMER.tar.gz",
    params:
        giggleexec = config["annotation"]["giggle"],
        gigglepath = config["annotation"]["giggledb"],
    threads:
	config["options"]["cores"]
    shell:
        "sort -k5nr {input.peak} | head -500 > {output.peaktop500};"
        "findMotifsGenome.pl {output.peaktop500} hg38 {homeroutput} -size given -mask -p {threads};"
	"tar -zcf {homeroutputgz} {homeroutput}"
	