process SPECIES_PROFILING_ONT {
	label "kraken2"
	tag { "Autodetecting reference using pooled samples" }

	publishDir "$outdir", overwrite: true, mode: 'copy'
	cpus 4

	input:
	tuple val(sample), path(reads)
	val(outdir)
	val(results)

	output:
	tuple val(sample), path("${sample}.kraken.report"), emit: report

	script:
	"""
		READS_DIR=\$(dirname $sample)
		
		# Pool Samples
		zcat \$READS_DIR/*.fastq.gz > pooled_samples.fastq.gz

		kraken2 -db $params.db --report ${sample}.kraken.report \
			--output - \
			--threads $task.cpus --memory-mapping \
			pooled_samples.fastq.gz
	"""

	stub:
	"""
		touch ${sample}.kraken.report
	"""
}


process DOWNLOAD_FASTA {
	label "kraken2"
	tag { "Downloading reference genome " }

	publishDir "$outdir", overwrite: true, mode: 'copy'
	
	input:
	tuple val(filename), path(kraken2_report)
	val(outdir)

	output:
	tuple val("reference"), path("reference.fa"), emit: genome
	tuple val("reference"), path("reference.gb"), emit: genbank
	tuple val("reference"), path("reference.log"), emit: log

	script:
	"""
		python3 $params.bin/download_reference.py $kraken2_report > reference.log

		# Unzip
		gunzip reference.fna.gz

		# Rename
		mv reference.fna reference.fa
	"""

	stub:
	"""
		touch reference.fa reference.gb reference.log
	"""
}


process EXTRACT_REFERENCE {
	label "kraken2"
	tag { "Extracting reference genome " }

	publishDir "$outdir", overwrite: true, mode: 'copy'
	
	input:
	path(filename)
	val(outdir)

	output:
	tuple val("reference"), path("reference.fa"), emit: genome

	script:
	"""
		python3 $params.bin/genbank2fasta.py $filename
		mv *.fa reference.fa
	"""

	stub:
	"""
		touch reference.fa
	"""
}
