configfile: "/home/u1357/encode/snakerbp/config_refv19_hepg2.yaml"
# configfile: "/home/u1357/encode/snakerbp/config_refv19_k562.yaml"

sampleFile = config["sampleFile"]
fileDir = config["fileDir"]
genomeFa = config["genomeFa"]
genomeGTF = config["genomeGTF"]
refDir = config["refDir"]


import pandas as pd
data = pd.read_table(sampleFile,header=0,sep="\t")
data.columns= ['experiment','type',"rep","sample","epName"]
SAMPLES = data.set_index("sample", drop=False)['sample'].to_dict()


rule all:
    input:
        (refDir + "star_index_1/sjdbList.fromGTF.out.tab"),
        expand(fileDir + "fastp/{sample}_fastp_end1.fastq.gz", sample=SAMPLES),
        expand(fileDir + "fastp/{sample}_fastp_end2.fastq.gz", sample=SAMPLES),
        expand(fileDir + "star_mapping_1/{sample}.SJ.out.tab", sample=SAMPLES),
        expand(fileDir + "star_mapping_2/{sample}.Aligned.out.bam", sample=SAMPLES),
        expand(fileDir + "sorted/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(fileDir + "rmats/{kdSample}/fromGTF.SE.txt", kdSample=kdSAMPLES),
        expand(fileDir + "rmatsDn/{kdSample}/sig.A3SS.MATS.JC.txt",kdSample=kdSAMPLES),
        expand(fileDir + "counts/{sample}.count", sample=SAMPLES),
        expand(fileDir + "merge/{kdSample}.counts.merge.txt", kdSample=kdSAMPLES,rep=REP),
        expand(fileDir + "deg/{kdSample}_diff_gene.txt", kdSample=kdSAMPLES),
        expand(fileDir + "bigwig/{sample}.Aligned.sortedByCoord.out.bw", sample=SAMPLES),
        
rule fastp:
    input:
        IP1 = fileDir + "raw/{sample}_end1.fastq.gz",
        IP2 = fileDir + "raw/{sample}_end2.fastq.gz",
    output:
        OP1 = fileDir + "fastp/{sample}_fastp_end1.fastq.gz",
        OP2 = fileDir + "fastp/{sample}_fastp_end2.fastq.gz",
        rpj = fileDir + "fastp/{sample}.json",
        rph = fileDir + "fastp/{sample}.html",
    log:
        fileDir + "log/fastp/{sample}.fastp.log",
    params:
        threads = 6,
    shell:
        """
        fastp -i {input.IP1} -o {output.OP1} -I {input.IP2} -O {output.OP2} -w {params.threads} -j {output.rpj} -h {output.rph} >{log} 2>&1
        rm {input.IP1} {input.IP2}
        """


# rule STAR_index_1:
#     input:
#         IP1 = genomeFa,
#         IP2 = genomeGTF
#     output:
#         LK = refDir + "star_index_1/sjdbList.fromGTF.out.tab"
#     threads: 20
#     log:
#         refDir + "star_index_1.log",
#     params:
#         genomeDir = refDir + "star_index_1/"
#     shell:
#         "STAR \
#         --runMode genomeGenerate \
#         --runThreadN 20 \
#         --genomeDir {params.genomeDir} \
#         --sjdbOverhang 100 \
#         --genomeFastaFiles {input.IP1} \
#         --sjdbGTFfile {input.IP2} >{log} 2>&1"


rule STAR_mapping_1:
    input:
        IP1 = fileDir + "fastp/{sample}_fastp_end1.fastq.gz",
	    IP2 = fileDir + "fastp/{sample}_fastp_end2.fastq.gz",
        LK = refDir + "star_index_1/sjdbList.fromGTF.out.tab",
    output:
        sjdb = fileDir + "star_mapping_1/{sample}.SJ.out.tab",
    log:
        fileDir + "log/star_mapping_1/{sample}.log"
    params:
        genomeDir = refDir + "star_index_1/",
        prefix = fileDir + "star_mapping_1/{sample}.",
    shell:
        "STAR \
        --genomeDir {params.genomeDir} \
	      --readFilesIn {input.IP1} {input.IP2} \
        --readFilesCommand gunzip -c \
        --runThreadN 12 \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 30000000000 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --sjdbScore 1 \
        --sjdbOverhang 100 \
        --outFilterMismatchNoverReadLmax  0.04 \
        --outFilterMismatchNmax 999 \
        --outFilterMultimapNmax 20 \
        --outSAMstrandField intronMotif \
        --outSAMmode None \
	      --outSAMtype None \
        --outFileNamePrefix {params.prefix} >{log} 2>&1"


## merge star_mapping_1 SJs
# nohup Rscript /home/u1357/encode/sjmerge/sjmerge.R /home/u1357/encode/rnaseq/hepg2/star_mapping_1/hepg2 >hepg2.log &
# nohup Rscript /home/u1357/encode/sjmerge/sjmerge.R /home/u1357/encode/rnaseq/k562/star_mapping_1/k562 >k562.log &


rule STAR_mapping_2:
    input:
        IP1 = fileDir + "fastp/{sample}_fastp_end1.fastq.gz",
        IP2 = fileDir + "fastp/{sample}_fastp_end2.fastq.gz",
        LK = refDir + "star_index_2/sjdbList.fromGTF.out.tab",
    output:
        bam = fileDir + "star_mapping_2/{sample}.Aligned.out.bam",
    log:
        fileDir + "log/star_mapping_2/{sample}.log",
    params:
        genomeDir = refDir + "star_index_2/",
        prefix = fileDir + "star_mapping_2/{sample}.",
    shell:
        """
        STAR \
        --genomeDir {params.genomeDir} \
	      --readFilesIn {input.IP1} {input.IP2} \
        --readFilesCommand gunzip -c \
        --runThreadN 8 \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 35000000000 \
        --alignEndsType EndToEnd \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --sjdbScore 1 \
        --sjdbOverhang 100 \
        --outFilterMismatchNoverReadLmax  0.04 \
        --outFilterMismatchNmax 999 \
        --outFilterMultimapNmax 20 \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS NM MD \
        --outSAMunmapped Within \
	      --outSAMheaderHD @HD VN:1.4 \
	      --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.prefix} >{log} 2>&1
        rm {input.IP1} {input.IP2}
        """


rule sort_index:
    input:
        fileDir + "star_mapping_2/{sample}.Aligned.out.bam"
    output:
        sort = fileDir + "sorted/{sample}.Aligned.sortedByCoord.out.bam",
        index = fileDir + "sorted/{sample}.Aligned.sortedByCoord.out.bam.bai",
    shell:
        """
        samtools sort -l 9 -@ 8 -o {output.sort} {input} 
        samtools index {output.sort} {output.index}
        rm {input}
        """


rule rMATS:
    input:
        nc = fileDir + "files/rmats/{kdSample}_nc.txt",
        kd = fileDir + "files/rmats/{kdSample}_kd.txt",
    output:
        LK = fileDir + "rmats/{kdSample}/fromGTF.SE.txt"
    log:
        fileDir + "log/rmats/{kdSample}.rmats.log",
    threads: 4
    params:
        outDir = fileDir + "rmats/{kdSample}/",
        tmp =  fileDir + "rmats/{kdSample}/tmp/",
        strand = "fr-firststrand",
        threads = 4,
    shell:
        "rmats.py --b1 {input.nc} --b2 {input.kd} --gtf {genomeGTF} -t paired --libType {params.strand} \
        --novelSS --readLength 100 --variable-read-length --nthread {params.threads} --od {params.outDir} --tmp {params.tmp} >{log} 2>&1"


rule rmatsDN:
    input:
        sample = fileDir + "rmats/{kdSample}/A3SS.MATS.JC.txt"
    output:
        fileDir + "rmatsDn/{kdSample}/sig.A3SS.MATS.JC.txt",
        fileDir + "rmatsDn/{kdSample}/total_ase_sig.txt",
    params:
        fileDir = fileDir + "rmats/{kdSample}/",
        outDir = fileDir + "rmatsDn/{kdSample}/",
        level = "nc,kd",
        nRep = 2,
        cutoff = 10,
        novel = "T",
    shell:
        "Rscript {codeDir}rmatsDN.R --fileDir {params.fileDir} --nRep {params.nRep} --cutoff {params.cutoff} --outDir {params.outDir} \
         --level {params.level} --novel {params.novel}"


rule featureCounts:
    input:
        bam = fileDir + "sorted/{sample}.Aligned.sortedByCoord.out.bam",
        index = fileDir + "sorted/{sample}.Aligned.sortedByCoord.out.bam.bai",
        gtf = genomeGTF,
    output:
        fileDir + "counts/{sample}.count",
    params:
        prefix = fileDir + "counts/{sample}",
        threads = 8,
        strand= 2, # 0
    shell:
        "Rscript script/run-featurecounts.R {input.bam} {input.gtf} {params.threads} {params.prefix} {params.strand}"


rule merge:
    input:
        IP1 = fileDir + "counts/{kdSample}_kd_rep1.count",
        IP2 = fileDir + "counts/{kdSample}_kd_rep2.count",
        IP3 = fileDir + "counts/{kdSample}_nc_rep1.count",
        IP4 = fileDir + "counts/{kdSample}_nc_rep2.count",
    output:
        fileDir + "merge/{kdSample}.counts.merge.txt",
    params:
        fileDir = fileDir + "counts/",
        kdsample = "{kdSample}",
        outDir = fileDir + "merge/",
    shell:
        "Rscript script/run-merge.R {params.fileDir} {params.kdsample} {params.outDir}"


rule DEG:
    input:
        counts = fileDir + "merge/{kdSample}.counts.merge.txt",
    output:
        fileDir + "deg/{kdSample}_diff_gene.txt",
    log:
        fileDir + "log/deg/{kdSample}.log",
    params:
        outDir = fileDir + "deg/{kdSample}",
    shell:
        "Rscript script/run-degEdgeR.R -c {input.counts} -g {gtfGene} --l1 {level1} --l2 {level2} -o {params.outDir} >{log} 2>&1"

rule bam2bw:
    input:
        bam = fileDir + "sorted/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bw = fileDir + "bigwig/{sample}.Aligned.sortedByCoord.out.bw",
    threads: 2
    shell:
        "bamCoverage -b {input.bam} -o {output.bw} -p 2"
