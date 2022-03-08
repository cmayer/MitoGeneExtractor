#Snakefile
configfile: "config.yaml"
rule all:
    input:
        expand("/home/usr/{sample}_out_alignment.fas", sample=config["samples"]),
        expand("/home/usr/{sample}_out_consensus.fas", sample=config["samples"])

rule fastq_dump_concat:
    input:
        "{sample}.sra"
    output:
        "{sample}_concat.fastq"
    shell:
        "module load sratoolkit/2.10.8 \n fastq-dump --split-e --readids {input} \n cat {wildcards.sample}*.fastq > {output}"

rule TrimGalore:
    input:
        "{sample}_concat.fastq"
    output:
        "{sample}_concat_trimmed.fq"
    shell:
        "perl /home/usr/bin/TrimGalore-0.6.6/trim_galore --no_report_file --dont_gzip --output_dir ./ {input}"

rule fasta2fas:
    input:
        "{sample}_concat_trimmed.fq"
    output:
        "{sample}_trimmed.fas"
    shell:
        "awk '(NR-1)%4 == 0 || (NR-2)%4==0' {input} | tr '@' '>' > {output}"

rule f2f_newname:
    input:
        "{sample}_trimmed.fas"
    output:
        "{sample}_trimmed_newname.fas"
    shell:
        "~/bin/fasta2fasta-v1.4 -r ' _' -u {input} {output}"

rule MitoGeneExtractor:
    input:
        DNA = "{sample}_trimmed_newname.fas",
        AA = "protein_consensus.fas"
    output:
        a = "{sample}_out_alignment.fas",
        b = "{sample}_out_consensus.fas",
        c = "{sample}_vulgar.txt",
        d = "{sample}_vulgar.txt.log"
    shell:
        "/home/usr/bin/MitoGeneExtractor -d {input.DNA} -p {input.AA} -o {output.a} -c {output.b} -V {output.c} -n 0 -t 0.5 -r 1 -e /home/usr/bin/exonerate"
