process CALCULATE_COVERAGE_DEPTH {
	label "genome_depth"
	publishDir "$outdir", overwrite: true, mode: 'copy'

	input:
	tuple val(barcode), path("${barcode}.mapped.bam")
	val(min_mq)
	val(outdir)

	output:
	tuple val(barcode), env(COV), emit: cov
	tuple val(barcode), env(DEPTH), emit: depth
	tuple val(barcode), path("${barcode}_depth.csv"), emit: summary
	tuple val(barcode), path("${barcode}_depth.tsv"), emit: tsv

	script:
	"""
	samtools depth -a --min-MQ $min_mq ${barcode}.mapped.bam > ${barcode}_depth.tsv
	python3 $params.bin/coverage_stats.py ${barcode}_depth.tsv ${barcode}
	COV=\$(tail -q -n1 coverage_stats.csv | cut -d ',' -f 10)
	DEPTH=\$(tail -q -n1 coverage_stats.csv | cut -d ',' -f 5)

	# Round to 0 barcodes with COV << 1 to avoid
	# inclusion, since BASH cannot handle them
	if [[ \$COV == *"e-"* ]]; then
		COV=0.0
		DEPTH=0.0
	fi

	mv coverage_stats.csv ${barcode}_depth.csv
	"""

	stub:
	"""
	touch ${barcode}_depth.tsv
	touch ${barcode}_depth.csv
	COV=0
	DEPTH=0
	MIN_READ_N=0
	"""
}


process RETRIEVE_CONTEXT {
	label "gen_consensus"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(variants), path(asmbl), path(reference)
	path(genbank)
	val(outdir)

	output:
	path("${barcode}_variants_context.tsv")

	script:
	"""
	python3 "$params.bin/snp_context_bcf.py" $variants $reference $genbank
	"""
	
	stub:
	"""
	touch ${barcode}_variants_context.tsv
	"""
}


process QC_ONT {
	label "quality_control"
	publishDir "$outdir", mode: 'copy'
	cpus 6

	input:
	tuple val(barcode), path(raw_reads)
	val(outdir)

	output:
	path("${barcode}_nanostats.txt"), emit: qc

	script:
	"""
	NanoPlot -t $task.cpus --N50 --fastq $raw_reads
	cp NanoStats.txt ${barcode}_nanostats.txt

	"""
	
	stub:
	"""
	mkdir figures
	touch ${barcode}_nanostats.txt 
	"""
}


process FILTER_READS {
	label "filter_reads"
	publishDir "$outdir", overwrite: true, mode: 'copy'
	cpus 6

	input:
	tuple val(barcode), path(raw_reads)
	val(Qthres)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}_filtered.fastq.gz"), emit: reads

	script:
	"""
	chopper --quality $Qthres --threads $task.cpus -i $raw_reads | gzip > ${barcode}_filtered.fastq.gz
	"""
	
	stub:
	"""
	touch ${barcode}_filtered.fastq.gz
	"""
}


process REFERENCE_STATS {
	label "gen_consensus"
	publishDir "$outdir", mode: 'copy'

	input:
	path(gb_assembly)
	val(outdir)

	output:
	tuple val("reference"), path("reference.stats"), emit: stats

	script:
	"""
	python3 "$params.bin/get_intergenic_space.py" $gb_assembly
	"""

	stub:
	"""
	touch reference.stats
	"""
}


process FIND_MLST {
	label "mlst"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(dn_assembly)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}_ST.tsv"), emit: mlst

	script:
	"""
	mlst --mincov 80 $dn_assembly > ${barcode}_ST.tsv
	"""

	stub:
	"""
	touch ${barcode}_ST.tsv
	"""
}


process FIND_wgMLST {
	label "mlst"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(dn_assembly)  // Kept to control timing
	val(denovo_path)
	val(reference)
	val(outdir)

	output:
	path("wgMLST.aln"), emit: msa
	path("wgMLST.aln"), emit: msa_wo_controls
	path("wgMLST_groups.tsv"), emit: groups

	script:
	"""
	# Import wgMLST database as wgmlst (default)
	if [ ! -f wgMLST.db ]; then
		wgmlst import wgmlst 'Staphylococcus aureus'
		mv wgmlst wgMLST.db
		cp wgMLST.db wgMLST_wo_controls.db
	fi

	# Add reference to the database
	wgmlst add --strain reference wgMLST.db $reference

	# Set limit on (approx) 85% size of reference to 125%
	REF_MINSIZE=\$(( \$(stat -c '%s' ${reference}) / 7 * 6 ))
	REF_MAXSIZE=\$(( \$(stat -c '%s' ${reference}) / 4 * 5 ))

	# Add samples to the database
	for ASMBL in $denovo_path/*.fa; do
		echo \$REF_MINSIZE \$REF_MAXSIZE
		if [ `stat --dereference -c '%s' \$ASMBL` -gt \$REF_MINSIZE ] && [ `stat --dereference -c '%s' \$ASMBL` -lt \$REF_MAXSIZE ]; then
			BARCODE=\$(basename \$ASMBL | cut -d "." -f 1)
			wgmlst add --strain \$BARCODE wgMLST.db \$ASMBL
			wgmlst add --strain \$BARCODE wgMLST_wo_controls.db \$ASMBL
		fi 
	done

	# Export alignment (just in case)
	wgmlst msa --output wgMLST.aln wgMLST.db
	wgmlst msa --output wgMLST_wo_controls.aln wgMLST_wo_controls.db

	# Export distance matrix - DO NOT! DISTANCE IN LOCI DIFFERENCES NOT SNP
	wgmlst distance -o wgMLST_distances.tsv wgMLST.db

	# Export clustering based on wgMLST
	wgmlst subgraph --export list --output wgMLST_groups.tsv wgMLST_distances.tsv
	"""

	stub:
	"""
	touch wgMLST.aln
	touch wgMLST_wo_controls.aln
	touch wgMLST_groups.tsv
	"""
}


process BED_FROM_BLAST {
	label "kraken2"
	publishDir "$outdir", mode: 'copy'
	
	input:
	val(reference)
	val(outdir)

	output:
	path("reference.bed"), emit: file
	path("repetitions.stats"), emit: stats

	script:
	"""
	python3 $params.bin/remove_repeated_regions.py $reference
	mv reference_db.stats repetitions.stats
	"""

	stub:
	"""
	touch reference.bed repetitions.stats
	"""
}
