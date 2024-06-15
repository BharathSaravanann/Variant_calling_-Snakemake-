rule all:
	input:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.vcf"


rule bwa_map:
    input:
        ref_file="/home/bharath/Lifecell/ref/hg38.fa",
        fastq1="/home/bharath/Lifecell/fastq/father_R1.fq.gz",
        fastq2="/home/bharath/Lifecell/fastq/father_R2.fq.gz"
    output:
        "output/bwa.sam"
    params:
        sample_id="father",
        rg_id="father",
        platform="ILLUMINA"
    shell:
        """
        bwa mem -t 4 -R "@RG\\tID:{params.rg_id}\\tPL:{params.platform}\\tSM:{params.sample_id}" {input.ref_file} {input.fastq1} {input.fastq2} > {output}
        """
	   
	   
rule sortbam:
	input:
	    "output/bwa.sam"
	    
	    
	output:
	    "output/file.sorted.bam"
	
	
	
	shell:
	    "gatk SortSam -I {input} -O {output} --SORT_ORDER coordinate"
	    
	    
rule indexbam:
	input:
	    "output/file.sorted.bam"
	
	
	
	output:
	    "output/file.sorted.idx"
	
	
	shell:
	    "gatk BuildBamIndex -I {input}"
	    
rule markduplicates:
	input:
	    "output/file.sorted.bam"
	
	
	output:
	    "output/file.sorted.dedup.bam"
	
	
	
	shell:
	    "gatk MarkDuplicatesSpark \
            -I {input} \
            -O {output}"
            
rule BQSR:
	input:
	    bamfile="output/file.sorted.dedup.bam",
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    knownsites="/home/bharath/Lifecell/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
	
	
	output:
	    "output/file.sorted.dedup.recal.data.table"
	
	
	
	shell:
	    "gatk BaseRecalibrator -I {input.bamfile} -R {input.ref_file} --known-sites {input.knownsites} -O {output}"
	  
	  
	  
rule applyBQSR:
	input:
	    bamfile="output/file.sorted.dedup.bam",
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    recal_data="output/file.sorted.dedup.recal.data.table"
	    
	
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.bam"
	
	
	shell:
	    "gatk ApplyBQSR -I {input.bamfile} -R {input.ref_file} --bqsr-recal-file {input.recal_data} -O {output}"
	    
	    
rule Haplotypecaller:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    bamfile="output/file.sorted.dedup.recal.bqsr.bam"
	    
	    
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.vcf"
	
	
	shell:
	    "gatk HaplotypeCaller -R {input.ref_file} -I {input.bamfile} -O {output}"
	

	    

