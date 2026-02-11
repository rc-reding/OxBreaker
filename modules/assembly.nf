process MAP_REFERENCE {
	label "reference_mapping"
	cpus 8

	publishDir "$outdir/", mode: 'copy'

	input:
	tuple val(barcode), path(reads_merged)
	path(ref_genome)
	path(ref_bed)
	val(min_mq)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.mapped.bam"), emit: bam

	script:
	"""
		# Index reference
		samtools faidx $ref_genome

		# Map reads to reference and filter not primary (256) unmapped (4) and supplt aligment reads
		minimap2 -x map-ont --seed $params.seed -N 0 -2 -Q -a --secondary=no -w 250 \
			 --sam-hit-only -t $task.cpus -o aln.sam $ref_genome $reads_merged;
		samtools view -F 2308 --min-MQ $min_mq --threads $task.cpus -b -o aln_filtered.bam aln.sam;
		samtools sort --threads $task.cpus -o ${barcode}.mapped.bam aln_filtered.bam;
	"""

	stub:
	"""
		touch ${barcode}.mapped.bam
	"""

}


process VARIANT_CALL_CLAIR3 {
	label "variant_call"
	cpus 8

	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(asmbl), env(asmbl_depth)
	val(min_freq)
	val(min_readN)
	val(min_mq)
	path(ref_genome)
	path(bed_file)
	val(outdir)
	
	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: vcf
	tuple val(barcode), path("${barcode}_report.vcf.gz"), emit: vcf_report

	script:
	"""
		samtools index $asmbl
		samtools faidx $ref_genome
		
		# Run clair3
		# Options min_coverage, min_mq, call_snp_only,
		# no_phasing_for_fa, are all tagged EXPERIMENTAL
		run_clair3.sh \
			--bam_fn=$asmbl \
			--ref_fn=$ref_genome \
			--bed_fn=$bed_file \
			--threads=$task.cpus \
			--platform="ont" \
			--include_all_ctgs \
			--no_phasing_for_fa \
			--call_snp_only \
			--min_coverage=$min_readN \
			--snp_min_af=$min_freq \
			--min_mq=$min_mq \
			--output=clair_out \
			--model_path=$params.model

		# Filter out INDELs ('--call_snp_only' is classed as EXPERIMENTAL in Clair3)
		bcftools view --exclude-types indels -Oz -o merge_output.vcf.gz clair_out/merge_output.vcf.gz
		tabix -p vcf merge_output.vcf.gz

		# Post-process (resolve gaps and null calls based on depth/freq data, and ambiguous multiallelic sites)
		python3 $params.bin/filterSNPs_clair3.py merge_output.vcf.gz $bed_file $min_freq $min_readN
		
		# Ensure headers consistent with vcf.gz files
		bcftools view -Oz -o ${barcode}.vcf.gz clair_out/merge_output_filtered.vcf.gz
		bcftools view -Oz -o ${barcode}_report.vcf.gz clair_out/merge_output_filtered_report.vcf.gz

	"""

	stub:
	"""
		touch ${barcode}.vcf.gz ${barcode}_report.vcf.gz
	"""
}


process VARIANT_CALL_BCFTOOLS {
	label "gen_consensus"
	cpus 8

	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(asmbl), env(asmbl_depth)
	val(min_freq)
	val(min_readN)
	val(min_mq)
	path(ref_genome)
	path(bed_file)
	val(outdir)
	
	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: vcf
	tuple val(barcode), path("${barcode}_report.vcf.gz"), emit: vcf_report

	script:
	"""
		samtools index $asmbl
		samtools faidx $ref_genome

		# Run bcftools pileup
		bcftools mpileup --threads $task.cpus \
			 --seed $params.seed \
			 --min-MQ $min_mq \
			 --fasta-ref $ref_genome \
			 --output-type b \
			 --write-index \
			 --regions-file $bed_file \
			 --output mpileup.vcf.gz $asmbl

		# Filter out INDELs (C > Python)
		bcftools view --exclude-types indels -Oz -o mpileup_filtered.vcf.gz mpileup.vcf.gz
		tabix -p vcf mpileup_filtered.vcf.gz

		# Run bcftools call
		bcftools call --threads $task.cpus \
			      --regions-file $bed_file \
			      --multiallelic-caller \
                              --ploidy 1 \
			      --output-type b \
			      --output variants.vcf.gz mpileup_filtered.vcf.gz

		# Post-process (resolve gaps and null calls based on depth/freq data, and ambiguous multiallelic sites)
		python3 $params.bin/filterSNPs_bcf.py variants.vcf.gz $bed_file $min_freq $min_readN
		mv variants_filtered.vcf.gz ${barcode}.vcf.gz
		mv variants_filtered_report.vcf.gz ${barcode}_report.vcf.gz
	"""

	stub:
	"""
		touch ${barcode}.vcf.gz ${barcode}_report.vcf.gz
	"""
}
