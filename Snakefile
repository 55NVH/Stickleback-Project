import os
import platform
import pandas as pd
import numpy as np
import string
import re, os
from itertools import chain

# # Import config file & parameters
# configfile: 'config.yaml'
def is_false(x):
    return x is None or x is False or (isinstance(x, str) and x.strip().lower() in ("false", "0", "none", ""))
# Import paths from config file
PATH = "{}".format(config['global']['snakemake'])
RESULTS_PATH = "{}".format(config['global']['results'])
GENOME_PATH = "{}/{}".format(config['genome']['path'], config['genome']['name'])
ANNOTATION_PATH = "{}/{}".format(config['annotation']['path'], config['annotation']['name'])
RIBO_INDEX_PATH = "{}/{}".format(config['alignment']['path'], config['alignment']['ribo'])
STAR_INDEX_PATH = "{}/{}".format(config['alignment']['path'], config['alignment']['star'])
SAMPLE_NAMES = list(config["samples"].keys())
FEATURECOUNTS_DIR = RESULTS_PATH + "/featurecounts"
MULTIQC_DIR = RESULTS_PATH + "/multiqc"
MULTIQC_HTML = MULTIQC_DIR + "/multiqc_report.html"
STAR_BAMS_DONE = RESULTS_PATH + "/done/star_bams_exist.ok"

DO_FASTQC_POST = not is_false(config.get("qc", {}).get("fastqc_post", True))
DO_BAMCOV = not is_false(config.get("qc", {}).get("bamcoverage", False))  # optional
DO_RIBO = not is_false(config.get("alignment", {}).get("ribo"))

def platform_info():
	import platform
	if platform.system() in "Darwin":
		return(".".join(("macOS", platform.machine())))
	else:
		return(".".join((platform.system().lower(), platform.machine())))


def strip_fastq_ext(fn):
    return re.sub(r"(\.fastq|\.fq)(\.gz)?$", "", fn)

def rename_out(inlist, suffix, outdir, outlist):
    if len(inlist) != len(outlist):
        raise ValueError(f"inlist/outlist length mismatch: {len(inlist)} vs {len(outlist)}")

    for infile, outfile in zip(inlist, outlist):
        base = os.path.basename(infile)
        stem = strip_fastq_ext(base)
        temp_out = os.path.join(outdir, stem + suffix)
        shell("mv {temp_out} {outfile}")

		
def expected_star_bams():
    """Final BAM+BAI that featureCounts should consume."""
    bams = []
    for sample in config["samples"]:
        ext = "PE" if len(config["samples"][sample]) == 2 else "SE"
        bams.append(f"{RESULTS_PATH}/mapped/{sample}_{ext}.Aligned.sortedByCoord.bam")
        bams.append(f"{RESULTS_PATH}/mapped/{sample}_{ext}.Aligned.sortedByCoord.bam.bai")
        # strand check file should also exist
        bams.append(f"{RESULTS_PATH}/mapped/{sample}_{ext}.StrandCheck.out.tab")
        # STAR final log (for MultiQC)
        bams.append(f"{RESULTS_PATH}/mapped/{sample}_{ext}.Log.final.out")
    return bams

def fc_strand_from_file(std_path):
    """For PE: convert StrandCheck to featureCounts -s value."""
    with open(std_path, "r") as f:
        val = f.readline().strip()
    if val == "FirstStrand":
        return "2"
    elif val == "SecondStrand":
        return "1"
    return "0"

def fastqc_post_targets():
    outs = []
    for sample in config["samples"]:
        if len(config["samples"][sample]) == 2:
            outs += [
                f"{RESULTS_PATH}/fastqc-post/{sample}_F_fastqc.html",
                f"{RESULTS_PATH}/fastqc-post/{sample}_R_fastqc.html",
            ]
        else:
            outs += [f"{RESULTS_PATH}/fastqc-post/{sample}_SE_fastqc.html"]
    return outs

rule all:
    input:
        RESULTS_PATH+"/sampleMetrics.tsv",
        RESULTS_PATH+"/sampleTableLong.tsv",
        # make fastqc-post actually run if enabled
        (fastqc_post_targets() if DO_FASTQC_POST else []),
        # featureCounts outputs for all samples
        expand(
            FEATURECOUNTS_DIR+"/{}_{}.featurecounts.txt".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        # optional bamcoverage outputs
        (expand(
            FEATURECOUNTS_DIR+"/{}_{}.bw".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ) if DO_BAMCOV else []),
        # MultiQC report
        MULTIQC_HTML

rule fastqc_pre_SE:
    input:
        fastq=lambda wc: config["samples"][wc.sample]
    output:
        html=RESULTS_PATH + "/fastqc-pre/{sample}_SE_fastqc.html",
        zip=temp(RESULTS_PATH + "/fastqc-pre/{sample}_SE_fastqc.zip")
    params:
        outdir=RESULTS_PATH + "/fastqc-pre"
    log:
        RESULTS_PATH + "/logs/fastqc-pre/{sample}.txt"
    threads: 1
    shell:
        r"""
        mkdir -p {params.outdir} {RESULTS_PATH}/logs/fastqc-pre
        fastqc -t {threads} --outdir {params.outdir} {input.fastq} &> {log}
        """

rule fastqc_pre_PE:
    input:
        fastq=lambda wc: config["samples"][wc.sample]   # [R1,R2]
    output:
        html1=RESULTS_PATH + "/fastqc-pre/{sample}_F_fastqc.html",
        html2=RESULTS_PATH + "/fastqc-pre/{sample}_R_fastqc.html",
        zip1=temp(RESULTS_PATH + "/fastqc-pre/{sample}_F_fastqc.zip"),
        zip2=temp(RESULTS_PATH + "/fastqc-pre/{sample}_R_fastqc.zip")
    params:
        outdir=RESULTS_PATH + "/fastqc-pre"
    log:
        RESULTS_PATH + "/logs/fastqc-pre/{sample}.txt"
    threads: 2
    shell:
        r"""
        mkdir -p {params.outdir} {RESULTS_PATH}/logs/fastqc-pre
        fastqc -t {threads} --outdir {params.outdir} {input.fastq} &> {log}
        """

rule trimgalore_SE:
	input:
		# html=expand(RESULTS_PATH+"/fastqc-pre/{{sample}}_{ext}_fastqc.html", ext=["SE"]),
		fastq=lambda wildcards: config["samples"][wildcards.sample]
	output:
		fastq=temp(expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}_trimmed.fq.gz", ext=["SE"]))
	params:
		basic="-q 20 --gzip --length 16 --no_report_file",
		outdir=RESULTS_PATH+"/trim_galore",
		temp=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}.fq.gz", ext=["SE"]),
		cluster='-N 1 -c 1 --mem=16G -t 160:00:00 -o logs/trimgalore.%A.{sample}.log'
	log:
		RESULTS_PATH+"/logs/trim_galore/{sample}_SE.fq.gz_trimming_report.txt"
	threads:
		2
	run:
		if config['barcode']:
			barcode = config['barcode']
			scriptfile=PATH+"/scripts/UMItoReadID.py"
			shell("python3 {scriptfile} -i {input.fastq} -o {params.temp} -b {barcode} -v &> {log}")
			shell("trim_galore {params.basic} --cores {threads} -o {params.outdir} {params.temp} &>> {log}")
			shell("rm {params.temp}")
		else:
			shell("trim_galore {params.basic} --cores {threads} -o {params.outdir} {input.fastq} &> {log}")
			rename_out(input["fastq"], ["_trimmed.fq.gz"], params["outdir"], output['fastq'])
			# rename_out(input["fastq"], [".fq.gz_trimming_report.txt"], RESULTS_PATH+"/logs/trim_galore", log['log'])

rule trimgalore_PE:
	input:
		# html=expand(RESULTS_PATH+"/fastqc-pre/{{sample}}_{ext}_fastqc.html", ext=["F", "R"]),
		fastq=lambda wildcards: config["samples"][wildcards.sample]
	output:
		fastq=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}_trimmed.fq.gz", ext=["F", "R"])
	params:
		basic="-q 20 --gzip --length 16 --no_report_file --paired --trim1",
		outdir=RESULTS_PATH+"/trim_galore",
		temp=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}.fq.gz", ext=["F", "R"]),
		cluster='-N 1 -c 1 --mem=16G -t 160:00:00 -o logs/trimgalore.%A.{sample}.log'
	log:
		RESULTS_PATH+"/logs/trim_galore/{sample}_PE.fq.gz_trimming_report.txt"
	threads:
		2
	run:
		def val_out(fq_path, suffix):
			"""
			Given an input fastq filename (fq.gz or fastq.gz),
			return Trim Galore output name with suffix like '_val_1.fq.gz' or '_val_2.fq.gz'.
			"""
			base = os.path.basename(fq_path)
			if base.endswith(".fastq.gz"):
		    		base = base[:-9]
			elif base.endswith(".fq.gz"):
				base = base[:-6]
			else:
				raise ValueError(f"Unexpected FASTQ extension: {fq_path}")
			return os.path.join(params["outdir"], base + suffix)
	
    	# Decide which input fastqs are passed to Trim Galore
		if config.get("barcode"):
			barcode = config["barcode"]
			scriptfile = PATH + "/scripts/UMItoReadID.py"

        # params.temp is a list of two temp fastqs: [*_F.fq.gz, *_R.fq.gz]
	     		shell("python3 {scriptfile} -i {input.fastq} -o {params.temp} -b {barcode} -v &> {log}")
	       		shell("trim_galore {params.basic} --cores {threads} -o {params.outdir} {params.temp} &>> {log}")

        # Trim Galore outputs are based on *each mate's basename*
		    	tg_r1 = val_out(params["temp"][0], "_val_1.fq.gz")
	       		tg_r2 = val_out(params["temp"][1], "_val_2.fq.gz")

	       		shell("mv {tg_r1} {output.fastq[0]}")
	       		shell("mv {tg_r2} {output.fastq[1]}")

	       		shell("rm -f {params.temp}")

		else:
			shell("trim_galore {params.basic} --cores {threads} -o {params.outdir} {input.fastq} &> {log}")

        # Trim Galore outputs are based on *each mate's basename*
			tg_r1 = val_out(input.fastq[0], "_val_1.fq.gz")
			tg_r2 = val_out(input.fastq[1], "_val_2.fq.gz")

			shell("mv {tg_r1} {output.fastq[0]}")
			shell("mv {tg_r2} {output.fastq[1]}")

if DO_RIBO:
	rule gunzip_rRNA_tRNA:
		input:
			RIBO_INDEX_PATH
		output:
			temp(os.path.splitext(RIBO_INDEX_PATH)[0])
		params:
			cluster='-N 1 -c 1 --mem=8G -t 160:00:00 -o logs/gunzip_genome.%A.log'
		shell:
			"""
			gunzip -c {input} > {output}
			"""

	rule index_rRNA_tRNA:
		input:
			os.path.splitext(RIBO_INDEX_PATH)[0] if os.path.splitext(RIBO_INDEX_PATH)[1] == ".gz" else RIBO_INDEX_PATH
		output:
			expand(RIBO_INDEX_PATH+".{ext}.bt2", ext=['1', '2', '3', '4', 'rev.1', 'rev.2'])
		params:
			basename=RIBO_INDEX_PATH,
			cluster='-N 1 -c 8 --mem=30G -t 160:00:00 -o logs/index_rRNA_tRNA.%A.log'
		log:
			RESULTS_PATH+"/logs/index_rRNA_tRNA/index_rRNA_tRNA.log"
		threads:
			8
		shell:
			"bowtie2-build --threads {threads} {input} {params.basename} &> {log}"

	rule remove_rRNA_tRNA_SE:
		input:
			fastq=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}_trimmed.fq.gz", ext=["SE"]),
			index1=expand(RIBO_INDEX_PATH+".{count}.bt2", count=['1', '2', '3', '4']),
			index2=expand(RIBO_INDEX_PATH+".rev.{count}.bt2", count=['1', '2'])
		output:
			fastq=temp(expand(RESULTS_PATH+"/rRNA_tRNA_removed/{{sample}}_{ext}_rRNA_tRNA_removed.fq.gz", ext=["SE"])),
			tmp=temp(expand(RESULTS_PATH+"/rRNA_tRNA_removed/{{sample}}_{tmp}.gz", tmp=["temp"])),
			metrics=RESULTS_PATH+"/rRNA_tRNA_removed/{sample}_SE_rRNA_tRNA_removed.txt"
		params:
			basic="--sensitive-local",
			basename=RIBO_INDEX_PATH,
			cluster='-N 1 -c 8 --mem=32G -t 160:00:00 -o logs/remove_rRNA_tRNA.%A.{sample}.log'
		log:
			RESULTS_PATH+"/logs/rRNA_tRNA_removed/{sample}_SE_stderr.log"
		threads:
			8
		run:
			shell("bowtie2 {params.basic} --threads {threads} --un-gz {output.fastq} -x {params.basename} -U {input.fastq} 2> {output.metrics} | gzip -c > {output.tmp} 2> {log}")

	rule remove_rRNA_tRNA_PE:
		input:
			fastq=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}_trimmed.fq.gz", ext=["F", "R"]),
			index1=expand(RIBO_INDEX_PATH+".{count}.bt2", count=['1', '2', '3', '4']),
			index2=expand(RIBO_INDEX_PATH+".rev.{count}.bt2", count=['1', '2'])
		output:
			fastq=temp(expand(RESULTS_PATH+"/rRNA_tRNA_removed/{{sample}}_{ext}_rRNA_tRNA_removed.fq.gz", ext=["F", "R"])),
			tmp=temp(expand(RESULTS_PATH+"/rRNA_tRNA_removed/{{sample}}_{tmp}.gz", tmp=["temp"])),
			metrics=RESULTS_PATH+"/rRNA_tRNA_removed/{sample}_PE_rRNA_tRNA_removed.txt"
		params:
			basic="--sensitive-local",
			basename=RIBO_INDEX_PATH,
			outname=RESULTS_PATH+"/rRNA_tRNA_removed/{sample}.fq.gz",
			outdir=RESULTS_PATH+"/rRNA_tRNA_removed",
			cluster='-N 1 -c 8 --mem=48G -t 160:00:00 -o logs/remove_rRNA_tRNA.%A.{sample}.log'
		log:
			RESULTS_PATH+"/logs/rRNA_tRNA_removed/{sample}_PE_stderr.log"
		threads:
			8
		run:
			fq1 = input["fastq"][0] 
			fq2 = input["fastq"][1] 
			shell("bowtie2 {params.basic} -p {threads} --un-conc-gz {params.outname} -x {params.basename} -1 {fq1} -2 {fq2} 2> {output.metrics} | gzip -c > {output.tmp} 2> {log}")
			inputList = [params["outname"].replace(".fq.gz", s) for s in [".fq.1.gz", ".fq.2.gz"]]
			# shell("python {}/scripts/pair_fastq_fast.py -l {} -r {}".format(PATH, *inputList))
			rename_out(inputList, [".fq.1.gz", ".fq.2.gz"], params["outdir"], output['fastq'])
else:
    rule skip_ribo_filter_PE:
        input:
            fwd = os.path.join(RESULTS_PATH, "trim_galore/{sample}_F_trimmed.fq.gz"),
            rev = os.path.join(RESULTS_PATH, "trim_galore/{sample}_R_trimmed.fq.gz"),
        output:
            fwd = os.path.join(RESULTS_PATH, "rRNA_tRNA_removed/{sample}_F_rRNA_tRNA_removed.fq.gz"),
            rev = os.path.join(RESULTS_PATH, "rRNA_tRNA_removed/{sample}_R_rRNA_tRNA_removed.fq.gz"),
        shell:
            r"""
            mkdir -p $(dirname {output.fwd})
            ln -sf $(realpath {input.fwd}) {output.fwd}
            ln -sf $(realpath {input.rev}) {output.rev}
            """

    rule ribo_metrics_stub_PE:
        input:
            fwd = os.path.join(RESULTS_PATH, "rRNA_tRNA_removed/{sample}_F_rRNA_tRNA_removed.fq.gz"),
            rev = os.path.join(RESULTS_PATH, "rRNA_tRNA_removed/{sample}_R_rRNA_tRNA_removed.fq.gz"),
        output:
            os.path.join(RESULTS_PATH, "rRNA_tRNA_removed/{sample}_PE_rRNA_tRNA_removed.txt")
        run:
            with open(output[0], "w") as f:
                f.write("rRNA/tRNA removal skipped (alignment.ribo: false)\n")
                f.write("No Bowtie2 filtering performed; using trimmed reads directly.\n")


rule fastqc_post_SE:
	input:
		fastq=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}_trimmed.fq.gz", ext=["SE"])
	output:
		html=expand(RESULTS_PATH+"/fastqc-post/{{sample}}_{ext}_fastqc.html", ext=["SE"]),
		gzip=expand(RESULTS_PATH+"/fastqc-post/{{sample}}_{ext}_fastqc.zip", ext=["SE"])
	params:
		outdir=RESULTS_PATH+"/fastqc-post",
		cluster='-N 1 -c 1 --mem=8G -t 160:00:00 -o logs/fastqc_post.%A.{sample}.log'
	log:
		RESULTS_PATH+"/logs/fastqc-post/{sample}.txt"
	threads:
		1
	run:
		shell("fastqc -t {threads} --outdir {params.outdir} {input.fastq} 2>&1 > {log}")
		rename_out(input["fastq"], ["_fastqc.html"], params["outdir"], output['html'])
		rename_out(input["fastq"], ["_fastqc.zip"], params["outdir"], output['gzip'])

rule fastqc_post_PE:
    input:
        fwd=os.path.join(RESULTS_PATH, "trim_galore/{sample}_F_trimmed.fq.gz"),
        rev=os.path.join(RESULTS_PATH, "trim_galore/{sample}_R_trimmed.fq.gz"),
    output:
        html_f=os.path.join(RESULTS_PATH, "fastqc-post/{sample}_F_fastqc.html"),
        html_r=os.path.join(RESULTS_PATH, "fastqc-post/{sample}_R_fastqc.html"),
        zip_f=os.path.join(RESULTS_PATH, "fastqc-post/{sample}_F_fastqc.zip"),
        zip_r=os.path.join(RESULTS_PATH, "fastqc-post/{sample}_R_fastqc.zip"),
    params:
        outdir=os.path.join(RESULTS_PATH, "fastqc-post"),
        cluster='-N 1 -c 2 --mem=8G -t 160:00:00 -o logs/fastqc_post.%A.{sample}.log'
    log:
        os.path.join(RESULTS_PATH, "logs/fastqc-post/{sample}.txt")
    threads: 2
    run:
        os.makedirs(params.outdir, exist_ok=True)

        # run fastqc on BOTH mates
        shell("fastqc -t {threads} --outdir {params.outdir} {input.fwd} {input.rev} &> {log}")

        # FastQC output names are based on input basenames:
        # e.g. FG_7dpa_repB_F_trimmed_fastqc.html
        f_base = os.path.basename(input.fwd).rsplit(".fq.gz", 1)[0]
        r_base = os.path.basename(input.rev).rsplit(".fq.gz", 1)[0]

        src_html_f = os.path.join(params.outdir, f_base + "_fastqc.html")
        src_zip_f  = os.path.join(params.outdir, f_base + "_fastqc.zip")
        src_html_r = os.path.join(params.outdir, r_base + "_fastqc.html")
        src_zip_r  = os.path.join(params.outdir, r_base + "_fastqc.zip")

        shell("mv {src_html_f} {output.html_f}")
        shell("mv {src_zip_f}  {output.zip_f}")
        shell("mv {src_html_r} {output.html_r}")
        shell("mv {src_zip_r}  {output.zip_r}")

rule gunzip_annotation:
	input:
		ANNOTATION_PATH
	output:
		os.path.splitext(ANNOTATION_PATH)[0]
	params:
		cluster='-N 1 -c 1 --mem=8G -t 160:00:00 -o logs/gunzip_annotation.%A.log'
	shell:
		"gunzip {input}> {output}"

rule gunzip_genome:
	input:
		GENOME_PATH
	output:
		os.path.splitext(GENOME_PATH)[0]
	params:
		cluster='-N 1 -c 1 --mem=8G -t 160:00:00 -o logs/gunzip_genome.%A.log'
	shell:
		"gunzip {input}> {output}"

rule index_star:
	input:
		genome=os.path.splitext(GENOME_PATH)[0] if os.path.splitext(GENOME_PATH)[1] == ".gz" else GENOME_PATH,
		annotation=os.path.splitext(ANNOTATION_PATH)[0] if os.path.splitext(ANNOTATION_PATH)[1] == ".gz" else ANNOTATION_PATH
	output:
		chrmlen=STAR_INDEX_PATH+"/chrNameLength.txt",
		SA=STAR_INDEX_PATH+"/SAindex"
	params:
		outdir=STAR_INDEX_PATH,
		overhang="--sjdbOverhang 99",
		cluster='-N 1 -c 8 --mem=64G -t 160:00:00 -o logs/index_star.%A.log'
	threads:
		16
	run:
		shell("STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.outdir} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.annotation} {params.overhang}")

def define_strand(counts, outfile):
	df = pd.read_table(counts, sep="\t", index_col=0, names=["Unstranded", "FirstStrand", "SecondStrand"])
	df = df.filter(regex="^(?!N_)", axis=0, ).sum()

	with open(outfile, "w") as of:
		if df['FirstStrand']/(df['FirstStrand']+df['SecondStrand']) > 0.55 :
			of.write("FirstStrand")
		elif df['FirstStrand']/(df['FirstStrand']+df['SecondStrand']) < 0.45 :
			of.write("SecondStrand")
		else:
			of.write("Unstranded")
	return(0)

rule map_star_SE:
	input:
		# html=expand(RESULTS_PATH+"/fastqc-post/{{sample}}_{ext}_fastqc.html", ext=["SE"]),
		fastq=expand(RESULTS_PATH+"/rRNA_tRNA_removed/{{sample}}_{ext}_rRNA_tRNA_removed.fq.gz", ext=["SE"]),
		chrmlen=STAR_INDEX_PATH+"/chrNameLength.txt"
	output:
		tmp=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.sortedByCoord.out.bam", ext=["SE"])),
		tmp2=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.toTranscriptome.out.bam", ext=["SE"])),
		# bam=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.sortedByCoord.bam", ext=["SE"])),
		cnt=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.ReadsPerGene.out.tab", ext=["SE"]),
		jct=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.SJ.out.tab", ext=["SE"]),
		sge=temp(directory(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}._STARgenome", ext=["SE"]))),
		spa=temp(directory(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}._STARpass1", ext=["SE"]))),
		std=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.StrandCheck.out.tab", ext=["SE"]),
		# wig=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Signal.Unique.str{num}.out.wig", ext=["SE"], num=["1", "2"])),
		# mul=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Signal.UniqueMultiple.str{num}.out.wig", ext=["SE"], num=["1", "2"])),
		log=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Log.final.out", ext=["SE"])
	params:
		genome=STAR_INDEX_PATH,
		outdir=RESULTS_PATH+"/mapped",
		alignment="--readFilesCommand zcat --genomeLoad NoSharedMemory --twopassMode Basic --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --limitBAMsortRAM 60000000000",
		rna="--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --quantMode TranscriptomeSAM GeneCounts",
		dna="--alignIntronMax 1 --alignMatesGapMax 300 --quantMode TranscriptomeSAM GeneCounts",
		# unstranded="--outSAMstrandField intronMotif",
		output="--outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterType BySJout --outSAMattributes All --outSAMtype BAM SortedByCoordinate",
		wiggle="--outWigType wiggle --outWigStrand Stranded --outWigNorm RPM",
		cluster='-N 1 -c 8 --mem=64G -t 160:00:00 -o logs/map_star.%A.{sample}.log'
	log:
		expand(RESULTS_PATH+"/logs/star_mapped/{{sample}}_{ext}.txt", ext=["SE"])
	threads:
		16
	run:
		prefix = "/".join((params["outdir"], os.path.basename(output["tmp"][0]).partition(".")[0]+"."))
		if config['seq_type'] == "RNAseq":
			shell("STAR --runMode alignReads --runThreadN {threads} --genomeDir {params.genome} --readFilesIn {input.fastq} --outFileNamePrefix {prefix} {params.alignment} {params.rna} {params.output} 2>&1 > {log}")
			define_strand(output["cnt"][0], output["std"][0])
		elif config['seq_type'] == "NETseq":
			shell("STAR --runMode alignReads --runThreadN {threads} --genomeDir {params.genome} --readFilesIn {input.fastq} --outFileNamePrefix {prefix} {params.alignment} {params.dna} {params.output} 2>&1 > {log}")
			define_strand(output["cnt"][0], output["std"][0])
		if config['barcode']:
			barcode = config['barcode']
			scriptfile=PATH+"/scripts/bamIDtoUMItag.py"
			shell("python3 {scriptfile} -i {output.tmp} -b {barcode} 2>&1 >> {log}")

rule map_star_PE:
	input:
		# html=expand(RESULTS_PATH+"/fastqc-post/{{sample}}_{ext}_fastqc.html", ext=["F", "R"]),
		fastq=expand(RESULTS_PATH+"/trim_galore/{{sample}}_{ext}_trimmed.fq.gz", ext=["F", "R"]),
		index=STAR_INDEX_PATH+"/SAindex",
		chrmlen=STAR_INDEX_PATH+"/chrNameLength.txt"
	output:
		tmp=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.sortedByCoord.out.bam", ext=["PE"])),
		tmp2=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.toTranscriptome.out.bam", ext=["PE"])),
		# bam=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.sortedByCoord.bam", ext=["PE"])),
		cnt=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.ReadsPerGene.out.tab", ext=["PE"]),
		jct=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.SJ.out.tab", ext=["PE"]),
		sge=temp(directory(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}._STARgenome", ext=["PE"]))),
		spa=temp(directory(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}._STARpass1", ext=["PE"]))),
		std=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.StrandCheck.out.tab", ext=["PE"]),
		# wig=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Signal.Unique.str{num}.out.wig", ext=["PE"], num=["1", "2"])),
		# mul=temp(expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Signal.UniqueMultiple.str{num}.out.wig", ext=["PE"], num=["1", "2"])),
		log=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Log.final.out", ext=["PE"])
	params:
		genome=STAR_INDEX_PATH,
		outdir=RESULTS_PATH+"/mapped",
		alignment="--readFilesCommand zcat --genomeLoad NoSharedMemory --twopassMode Basic --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --limitBAMsortRAM 60000000000",
		rna="--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --quantMode TranscriptomeSAM GeneCounts",
		dna="--alignIntronMax 1 --alignMatesGapMax 300 --quantMode TranscriptomeSAM GeneCounts",
		# unstranded="--outSAMstrandField intronMotif",
		output="--outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterType BySJout --outSAMattributes All --outSAMtype BAM SortedByCoordinate",
		wiggle="--outWigType wiggle --outWigStrand Stranded --outWigNorm RPM",
		cluster='-N 1 -c 8 --mem=64G -t 160:00:00 -o logs/map_star.%A.{sample}.log',
		#startmp=lambda wildcards:f"--outTmpDir /tmp/STAR_{wildcards.sample}"
	log:
		expand(RESULTS_PATH+"/logs/star_mapped/{{sample}}_{ext}.txt", ext=["PE"])
	threads:
		16 # STAR fails when using more than 16 threads
	run:
		prefix = "/".join((params["outdir"], os.path.basename(output["tmp"][0]).partition(".")[0]+"."))
		if config['seq_type'] == "RNAseq":
			shell("STAR --runMode alignReads --runThreadN {threads} --genomeDir {params.genome} --readFilesIn {input.fastq} --outFileNamePrefix {prefix} {params.alignment} {params.rna} {params.output} > {log} 2>&1")
			define_strand(output["cnt"][0], output["std"][0])
		elif config['seq_type'] == "NETseq":
			shell("STAR --runMode alignReads --runThreadN {threads} --genomeDir {params.genome} --readFilesIn {input.fastq} --outFileNamePrefix {prefix} {params.alignment} {params.dna} {params.output} 2>&1 > {log}")
			define_strand(output["cnt"][0], output["std"][0])
		# if config['barcode']:
		# 	barcode = config['barcode']
		# 	scriptfile=PATH+"/scripts/bamIDtoUMItag.py"
		# 	shell("python3 {scriptfile} -i {output.tmp} -b {barcode} 2>&1 >> {log}")
		# 	shell("python3 {scriptfile} -i {output.tmp2} -b {barcode} 2>&1 >> {log}")

rule fix_barcode_PE:
	input:
		bam1=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.sortedByCoord.out.bam", ext=["PE"]),
		bam2=expand(RESULTS_PATH+"/mapped/{{sample}}_{ext}.Aligned.toTranscriptome.out.bam", ext=["PE"])
	output:
		bam1=temp(expand(RESULTS_PATH+"/barcode_fixed/{{sample}}_{ext}.Aligned.sortedByCoord.out.bam", ext=["PE"])),
		bam2=temp(expand(RESULTS_PATH+"/barcode_fixed/{{sample}}_{ext}.Aligned.toTranscriptome.out.bam", ext=["PE"]))
	params:
		barcode = config['barcode'],
		scriptfile=PATH+"/scripts/bamIDtoUMItag.py"
	log:
		expand(RESULTS_PATH+"/logs/barcode_fixed/{{sample}}_{ext}.txt", ext=["PE"])
	threads:
		2
	run:
		shell("python3 {params.scriptfile} -i {input.bam1} -o {output.bam1} -b {params.barcode} 2>&1 >> {log}")
		shell("python3 {params.scriptfile} -i {input.bam2} -o {output.bam2} -b {params.barcode} 2>&1 >> {log}")

rule flag_xs_SE:
    input:
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.sortedByCoord.out.bam",
        # keep if you still want to also tag transcriptome bam (optional)
        txs=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.toTranscriptome.out.bam"
    output:
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.sortedByCoord.bam",
        bai=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.sortedByCoord.bam.bai",
        # OPTIONAL: if you still need transcriptome BAM downstream, keep these:
        bam2=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.toTranscriptome.bam",
        bai2=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.toTranscriptome.bam.bai",
        tmp=temp(f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.toTranscriptome.tmp.bam")
    params:
        cluster='-N 1 -c 8 --mem=32G -t 160:00:00 -o logs/flag_xs.%A.{sample}.log'
    log:
        f"{RESULTS_PATH}/logs/flag_xs/{{sample}}_SE.txt"
    threads: 4
    run:
        os.makedirs(f"{RESULTS_PATH}/logs/flag_xs", exist_ok=True)
        scriptfile = PATH + "/scripts/tagXSstrandedData.awk"

        # coordinate bam
        shell(
            "samtools view --threads {threads} -h {input.bam} 2> {log} "
            "| awk -f {scriptfile} 2>> {log} "
            "| samtools view --threads {threads} -b - 2>> {log} > {output.bam}"
        )
        shell("samtools index {output.bam} 2>> {log}")

        # OPTIONAL transcriptome bam (remove this whole block if you don't need it)
        shell(
            "samtools view --threads {threads} -h {input.txs} 2>> {log} "
            "| awk -f {scriptfile} 2>> {log} "
            "| samtools view --threads {threads} -b - 2>> {log} > {output.tmp}"
        )
        shell("samtools sort --threads {threads} -O BAM -o {output.bam2} {output.tmp} 2>> {log}")
        shell("samtools index {output.bam2} 2>> {log}")

rule flag_xs_PE:
    input:
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_PE.Aligned.sortedByCoord.out.bam",
        std=f"{RESULTS_PATH}/mapped/{{sample}}_PE.StrandCheck.out.tab"
    output:
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_PE.Aligned.sortedByCoord.bam",
        bai=f"{RESULTS_PATH}/mapped/{{sample}}_PE.Aligned.sortedByCoord.bam.bai"
    params:
        cluster='-N 1 -c 8 --mem=32G -t 160:00:00 -o logs/flag_xs.%A.{sample}.log'
    log:
        f"{RESULTS_PATH}/logs/flag_xs/{{sample}}_PE.txt"
    threads: 4
    run:
        os.makedirs(f"{RESULTS_PATH}/logs/flag_xs", exist_ok=True)

        # StrandCheck decides how your awk tags XS (same logic as before)
        with open(input["std"], "r") as infile:
            first_line = infile.readline().strip()
            if first_line == "FirstStrand":
                strand = 1
            elif first_line == "SecondStrand":
                strand = 2
            else:
                strand = 0

        scriptfile = PATH + "/scripts/tagXSstrandedData.awk"

        shell(
            "samtools view --threads {threads} -h {input.bam} 2> {log} "
            "| awk -v strType={strand} -f {scriptfile} 2>> {log} "
            "| samtools view --threads {threads} -b - 2>> {log} > {output.bam}"
        )
        shell("samtools index {output.bam} 2>> {log}")


rule star_bams_exist_checkpoint:
    """
    Hard gate: ensure ALL final STAR BAM/BAI/StrandCheck/Log.final.out exist
    before any featureCounts (prevents partial runs / filesystem race).
    """
    input:
        expected_star_bams()
    output:
        touch(STAR_BAMS_DONE)
    run:
        missing = [p for p in input if not os.path.exists(p)]
        if missing:
            raise ValueError("Missing post-STAR files (cannot start featureCounts):\n" + "\n".join(missing))

# ------------------------------
# featureCounts: SE fixed -s 2
# ------------------------------
rule featurecounts_SE:
    input:
        gate=STAR_BAMS_DONE,
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_SE.Aligned.sortedByCoord.bam",
        gtf=os.path.splitext(ANNOTATION_PATH)[0] if os.path.splitext(ANNOTATION_PATH)[1] == ".gz" else ANNOTATION_PATH
    output:
        txt=f"{FEATURECOUNTS_DIR}/{{sample}}_SE.featurecounts.txt",
        summary=f"{FEATURECOUNTS_DIR}/{{sample}}_SE.featurecounts.txt.summary"
    params:
        outdir=FEATURECOUNTS_DIR,
        cluster='-N 1 -c 8 --mem=16G -t 160:00:00 -o logs/featurecounts.%A.{sample}.log'
    threads: 8
    log:
        f"{RESULTS_PATH}/logs/featurecounts/{{sample}}_SE.txt"
    run:
        shell(
            "featureCounts -T {threads} -s 2 -t exon -g gene_id "
            "-a {input.gtf} -o {output.txt} {input.bam} &> {log}"
        )

# ------------------------------
# featureCounts: PE strand-aware
# ------------------------------
rule featurecounts_PE:
    input:
        gate=STAR_BAMS_DONE,
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_PE.Aligned.sortedByCoord.bam",
        gtf=os.path.splitext(ANNOTATION_PATH)[0] if os.path.splitext(ANNOTATION_PATH)[1] == ".gz" else ANNOTATION_PATH,
        std=f"{RESULTS_PATH}/mapped/{{sample}}_PE.StrandCheck.out.tab"
    output:
        txt=f"{FEATURECOUNTS_DIR}/{{sample}}_PE.featurecounts.txt",
        summary=f"{FEATURECOUNTS_DIR}/{{sample}}_PE.featurecounts.txt.summary"
    params:
        outdir=FEATURECOUNTS_DIR,
        cluster='-N 1 -c 8 --mem=16G -t 160:00:00 -o logs/featurecounts.%A.{sample}.log'
    threads: 8
    log:
        f"{RESULTS_PATH}/logs/featurecounts/{{sample}}_PE.txt"
    run:
        os.makedirs(params.outdir, exist_ok=True)
        strand = fc_strand_from_file(input.std)
        shell(
            "featureCounts -T {threads} -s 2 -p -t exon -g gene_id "
            "-a {input.gtf} -o {output.txt} {input.bam} &> {log}"
        )

# ------------------------------
# Optional bamCoverage (like featurecount.sh)
# ------------------------------
rule bamcoverage_star_bam:
    input:
        gate=STAR_BAMS_DONE,
        bam=f"{RESULTS_PATH}/mapped/{{sample}}_{{ext}}.Aligned.sortedByCoord.bam"
    output:
        bw=f"{FEATURECOUNTS_DIR}/{{sample}}_{{ext}}.bw"
    threads: 8
    log:
        f"{RESULTS_PATH}/logs/bamcoverage/{{sample}}_{{ext}}.txt"
    shell:
        r"""
        mkdir -p {FEATURECOUNTS_DIR} {RESULTS_PATH}/logs/bamcoverage
        bamCoverage -b {input.bam} -o {output.bw} --binSize 10 --normalizeUsing RPKM -p {threads} &> {log}
        """

# ------------------------------
# MultiQC over trim/fastqc/star/featurecounts
# ------------------------------
rule multiqc:
    input:
        gate=STAR_BAMS_DONE
    output:
        MULTIQC_HTML
    params:
        outdir=MULTIQC_DIR
    log:
        RESULTS_PATH + "/logs/multiqc.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        conda run -n rnaseq_multiqc multiqc {RESULTS_PATH} -o {params.outdir} &> {log}
        """


rule metrics:
    input:
        trimmed=expand(
            RESULTS_PATH+"/logs/trim_galore/{}_{}.fq.gz_trimming_report.txt".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        mapped=expand(
            RESULTS_PATH+"/mapped/{}_{}.Log.final.out".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        strand=expand(
            RESULTS_PATH+"/mapped/{}_{}.StrandCheck.out.tab".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        fc_summary=expand(
            FEATURECOUNTS_DIR+"/{}_{}.featurecounts.txt.summary".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        )
    output:
        RESULTS_PATH+"/sampleMetrics.tsv"
    params:
        cluster='-N 1 -c 1 --mem=4GB -t 160:00:00 -o logs/report.%A.log'
    run:
        idsList = [s for s in config["samples"]]

        # Trim Galore: "sequences processed in total" and "shorter than the length cutoff"
        iniList = [x for x in shell(
            "for infile in {input.trimmed}; do grep -m 1 'sequences processed in total' $infile | awk '{{print $1}}'; done",
            iterable=True
        )]
        triList = [x for x in shell(
            "for infile in {input.trimmed}; do grep -m 1 'shorter than the length cutoff' $infile | awk '{{print $(NF-1)}}'; done",
            iterable=True
        )]

        # STAR mapping
        uniList = [x for x in shell(
            "for infile in {input.mapped}; do grep 'Uniquely mapped reads number' $infile | awk '{{print $NF}}'; done",
            iterable=True
        )]
        mulList = [x for x in shell(
            "for infile in {input.mapped}; do grep 'Number of reads mapped to too many loci' $infile | awk '{{print $NF}}'; done",
            iterable=True
        )]
        strList = [x for x in shell(
            "for infile in {input.strand}; do awk '{{print $1}}' $infile; done",
            iterable=True
        )]

        # featureCounts summary: Assigned
        assigned = [x for x in shell(
            "for infile in {input.fc_summary}; do grep -w '^Assigned' $infile | awk '{{print $2}}'; done",
            iterable=True
        )]

        DF = pd.DataFrame({
            "Total reads": iniList,
            "Trimmed (too short)": triList,
            "Uniquely mapped": uniList,
            "Multi-hits mapped": mulList,
            "Library strand": strList,
            "featureCounts Assigned": assigned
        }, index=idsList)

        DF.to_csv(output[0], sep="\t", index_label="Sample")

rule report:
    input:
        bam=expand(
            RESULTS_PATH+"/mapped/{}_{}.Aligned.sortedByCoord.bam.bai".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        strand=expand(
            RESULTS_PATH+"/mapped/{}_{}.StrandCheck.out.tab".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        fc=expand(
            FEATURECOUNTS_DIR+"/{}_{}.featurecounts.txt".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        )
    output:
        RESULTS_PATH+"/sampleTable.tsv"
    params:
        cluster='-N 1 -c 1 --mem=4GB -t 160:00:00 -o logs/report.%A.log'
    run:
        idsList = [s for s in config["samples"]]
        bamList = [os.path.basename(p).replace(".bai", "") for p in input["bam"]]
        stdList = [os.path.basename(p) for p in input["strand"]]
        fcList  = [os.path.basename(p) for p in input["fc"]]

        DF = pd.DataFrame({
            "BamFile": bamList,
            "GeneCount": fcList,
            "LibraryStrand": stdList
        }, index=idsList)

        DF.to_csv(output[0], sep="\t", index_label="Sample")

rule report_long:
    input:
        bam=expand(
            RESULTS_PATH+"/mapped/{}_{}.Aligned.sortedByCoord.bam.bai".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        strand=expand(
            RESULTS_PATH+"/mapped/{}_{}.StrandCheck.out.tab".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        fc=expand(
            FEATURECOUNTS_DIR+"/{}_{}.featurecounts.txt".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ),
        bw=(expand(
            FEATURECOUNTS_DIR+"/{}_{}.bw".format(sample, ext)
            for sample in config["samples"]
            for ext in (["PE"] if len(config["samples"][sample])==2 else ["SE"])
        ) if DO_BAMCOV else [])
    output:
        RESULTS_PATH+"/sampleTableLong.tsv"
    params:
        cluster='-N 1 -c 1 --mem=4GB -t 160:00:00 -o logs/report.%A.log'
    run:
        idsList = [s for s in config["samples"]]
        bamList = [os.path.basename(p).replace(".bai", "") for p in input["bam"]]
        stdList = [os.path.basename(p) for p in input["strand"]]
        fcList  = [os.path.basename(p) for p in input["fc"]]

        data = {
            "BamFile": bamList,
            "GeneCount": fcList,
            "LibraryStrand": stdList
        }

        if DO_BAMCOV:
            bwList = [os.path.basename(p) for p in input["bw"]]
            data["BigWig"] = bwList

        DF = pd.DataFrame(data, index=idsList)
        DF.to_csv(output[0], sep="\t", index_label="Sample")

