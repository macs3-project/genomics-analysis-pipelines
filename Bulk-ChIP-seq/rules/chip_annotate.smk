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
