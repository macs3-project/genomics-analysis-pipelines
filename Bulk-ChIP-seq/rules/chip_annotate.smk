# qc for sequencing alignment from each replicate
rule chip_peakanno:
    input:
        peak = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak",
    output:
        pdf = "{OUT_DIR}/Analysis/{name}_peaks.narrowPeak_chip_annotation.pdf"
    params:
        gmtfolder = config["annotation"]["GMT"],
    shell:
        "utils/chip_annotate.R {input.peak} {params.gmtfolder};"
