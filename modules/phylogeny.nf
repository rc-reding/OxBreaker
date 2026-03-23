process GENERATE_CONSENSUS {
	label "genome_alignments"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), \
	      path(barcode_variant), \
	      path(assembly), \
	      path(barcode_depth), \
	      env(asmbl_depth)
	path(ref_genome)
	path(bed_file)
	val(min_depth)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.fa"), emit: fa

	script:
	"""
	FNAME=\$(echo $barcode | cut -d'.' -f1)

	awk -v MIN_DEPTH=$min_depth '{if (\$3 < $min_depth) print \$1"\t"\$2}' $barcode_depth > low_coverage.tsv
	
	# Index vcf
	cp $barcode_variant tmp.vcf.gz
	gunzip tmp.vcf.gz
	bgzip -o ${barcode}.vcf.gz tmp.vcf
	tabix -p vcf ${barcode}.vcf.gz

	# Gsample - file
	bcftools consensus --samples - \
			   --mask low_coverage.tsv \
			   --mark-del - \
			   --mark-snv lc \
			   --fasta-ref $ref_genome \
			   --output ${barcode}.fa \
			   ${barcode}.vcf.gz

	# Change sample name in fasta file
	sed -i -e "s/>/>\$FNAME /" "${barcode}.fa"
	"""

	stub:
	"""
	touch ${barcode}.fa
	"""
}


process CONSENSUS_2_MSA{
	publishDir "$outdir/", mode: 'copy'

	input:
	val(consensus_timed)  // no touching, helps control timing
	path(consensus)
	path(ref_genome)
	val(outdir)

	output:
	path("msa.fa"), emit: alignment
	path("msa_stats.txt"), emit: alignment_stats
	path("msa_w_controls.fa"), emit: alignment_w_controls

	script:
	"""
		if [ ! -f msa_w_controls.fa ]; then
			cat $ref_genome | sed 's/>/>reference /' > msa_w_controls.fa
		fi

		# Avoids ERR if no sample eligible
		if [ ! -f msa.fa ]; then
			touch msa.fa
			echo "Sample Name" \\\t "Proportion of genome masked" > msa_stats.txt
		fi

		for SEQ in \$(ls $consensus/*.fa); do
			# If most of the assembly is not masked...
			Ns=\$(tail -n +2 \$SEQ | grep "N" | wc -c)
			ns=\$(tail -n +2 \$SEQ | grep "n" | wc -c)
			gaps=\$(tail -n +2 \$SEQ | grep "-" | wc -c)

			GENOME_SIZE=\$(tail -n +2 \$SEQ | wc -c)
			PROPORTION_MASKED=\$(echo "scale=2; (\$Ns + \$ns + \$gaps) / \$GENOME_SIZE" | bc)

			echo \$(basename \$SEQ) \\\t \$PROPORTION_MASKED >> msa_stats.txt 


			if [ \$(echo "\$PROPORTION_MASKED <= $params.masking_tolerance" | bc) -eq 1 ]; then
				cat \$SEQ >> msa_w_controls.fa
				cat \$SEQ >> msa.fa
			fi
		done

	"""

	stub:
	"""
		touch msa.fa msa_w_controls.fa msa_stats.txt
	"""
	
}


process PREALIGN_GENOMES {
	label "genome_alignments"
	publishDir "$outdir/", mode: 'copy'

	input:
	tuple val(barcode), env(READ_BREADTH)
	val(barcode_nil)  // DO NOT DELETE: USED TO ENSURE PREV STEP COMPLETED
	path(path_var)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.vcf.gz"), emit: preprocessed
	path("${barcode}.vcf.gz.csi"), emit: index

	script:
	MIN_BREADTH=0.7
	"""
	if [[ \$READ_BREADTH > $MIN_BREADTH || \$READ_BREADTH == $MIN_BREADTH ]]; then
		if [[ `stat -c '%s' $path_var/${barcode}.vcf.gz` -gt 0 ]]; then
			cp vcf/${barcode}.vcf.gz ${barcode}.vcf.gz

			gunzip ${barcode}.vcf.gz
			bgzip ${barcode}.vcf
			tabix -p vcf --csi ${barcode}.vcf.gz
		else
			touch ${barcode}.vcf.gz
			touch ${barcode}.vcf.gz.csi
		fi
	else
		touch ${barcode}.vcf.gz
		touch ${barcode}.vcf.gz.csi
	fi
	"""

	stub:
	"""
	touch ${barcode}.vcf.gz
	touch ${barcode}.vcf.gz.csi
	"""
}


process GET_SAMPLES_NUMBER {
	input:
	val(barcode_nil)  // DO NOT DELETE: USED TO ENSURE PREV STEP COMPLETED
	path(path_var)

	output:
	env(sampleN)

	script:
	"""
		sampleN=\$(ls $path_var/*.gz |  wc -l)
	"""

	stub:
	"""
		sampleN=100
	"""
}


process GENERATE_PHYLOGENY {
	label "phylogeny_gubbins"
	cpus 8
	publishDir "$outdir", mode: 'copy'

	input:
	path("msa.fa")
	val(outdir)

	output:
	path("recomb_free.aln"), emit: aln
	path("corrected_tree.nwk"), emit: tree

	script:
	"""
	run_gubbins.py --prefix gubbins --seed $params.seed --model GTRGAMMA --recon-model GTRGAMMA \
		--min-snps 3 --p-val 0.01 --extensive-search --min-window-size 100 \
		--sh-test --threads $task.cpus --iterations 10 msa.fa

	mask_gubbins_aln.py --aln gubbins.filtered_polymorphic_sites.fasta \
		--gff gubbins.recombination_predictions.gff \
		--out recomb_free.aln

	# mv gubbins.filtered_polymorphic_sites.fasta recomb_free.aln
	mv gubbins.final_SH_support_tree.tre corrected_tree.nwk
	"""

	stub:
	"""
	touch recomb_free.aln
	touch corrected_tree.nwk
	"""
}


process GENERATE_PHYLOGENY_RAxML {
	label "phylogeny"
	cpus 8
	publishDir "$outdir", mode: 'copy'

	input:
	path("msa.fa")
	val(outdir)

	output:
	path("sh_tree.nwk"), emit: tree

	script:
	"""
	# Look for best tree (no bootstraping)
	raxmlHPC-PTHREADS-SSE3 -T $task.cpus -m GTRGAMMA \
			  -o reference -s msa.fa \
			  -n BTREE -k -p 100100

	# Run SH-like test on best tree
	raxmlHPC-PTHREADS-SSE3 -T $task.cpus -m GTRGAMMA \
			  -o reference -s msa.fa \
			  -n FINAL -p 011011 \
			  -f J -t RAxML_bestTree.BTREE
	
	mv RAxML_fastTreeSH_Support.FINAL sh_tree.nwk
	"""

	stub:
	"""
	touch sh_tree.nwk
	"""
}


process CORRECT_RECOMBINATION {
	label "phylogeny"
	publishDir "$outdir", mode: 'copy'

	input:
	path(best_tree)
	path("msa.fa")
	val(outdir)

	output:
	path("recomb_free.aln"), emit: aln
	path("corrected_tree.nwk"), emit: tree

	script:
	"""
	ClonalFrameML $best_tree msa.fa corrected_tree
	mv corrected_tree.labelled_tree.newick  corrected_tree.nwk
	mv corrected_tree.ML_sequence.fasta recomb_free.aln
	"""

	stub:
	"""
	touch recomb_free.aln
	touch corrected_tree.nwk
	"""
}


process DISTANCE_MATRIX {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	path(msa)
	val(outdir)

	output:
	path("distance_matrix.tsv"), emit: snp

	script:
	"""
	snp-dists -q $msa > distance_matrix.tsv
	"""

	stub:
	"""
	touch distance_matrix.tsv
	"""
}
