process GENERATE_CONSENSUS {
	label "genome_alignments"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(barcode_variant), path(barcode_depth), env(asmbl_depth)
	path(assembly)
	path(ref_genome)
	val(outdir)

	output:
	tuple val(barcode), path("${barcode}.fa"), emit: fa

	script:
	DEPTH_THRESHOLD=0.5
	"""
	FNAME=\$(echo $barcode | cut -d'.' -f1)
	MIN_DEPTH=\$(bc -s <<< "\$asmbl_depth * $DEPTH_THRESHOLD")
	MIN_DEPTH=10

	# Generate mask based on depth	
	awk -v MIN_DEPTH=\$MIN_DEPTH '{if (\$3 < MIN_DEPTH) print \$1"\t"\$2}' $barcode_depth > mask.txt
	
	# Index vcf
	tabix -p vcf $barcode_variant
	# Normalise indels/snps
	#bcftools norm --fasta-ref $ref_genome $barcode_variant -Oz -o variants.norm.vcf.gz
	# Index vcf
	#tabix -p vcf variants.norm.vcf.gz

	# Generate fasta file
	bcftools consensus --mask mask.txt \
					   --mark-del - \
					   --mark-snv lc \
					   --fasta-ref $ref_genome \
					   --output ${barcode}.fa \
					   $barcode_variant

	# Change sample name in fasta file
	sed -i "s/>/>\$FNAME /" "${barcode}.fa"
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
	path("msa_w_controls.fa"), emit: alignment_w_controls

	script:
	"""
		if [ ! -f msa_w_controls.fa ]; then
			cat $ref_genome | sed 's/>/>reference /' > msa_w_controls.fa
		fi

		# Avoids ERR if no sample eligible
		if [ ! -f msa.fa ]; then
			touch msa.fa
		fi

		for SEQ in \$(ls $consensus/*.fa); do
			# If most of the assembly is not masked...
			Ns=\$(tail -n +2 \$SEQ | grep "N" | wc -c)
			GENOME_SIZE=\$(tail -n +2 \$SEQ | wc -c)
			PROPORTION_MASKED=\$(echo "scale=2; \$Ns / \$GENOME_SIZE" | bc)

			if [ \$(echo "\$PROPORTION_MASKED < .25" | bc) -eq 1 ]; then
				cat \$SEQ >> msa_w_controls.fa
				if [[ ! \$(basename \$SEQ) == *"01"* && ! \$(basename \$SEQ) == *"02"* &&
		      		      ! \$(basename \$SEQ) == *"49"* && ! \$(basename \$SEQ) == *"50"* &&
		      		      ! \$(basename \$SEQ) == *"65"* && ! \$(basename \$SEQ) == *"66"* ]]; then
					cat \$SEQ >> msa.fa
				fi
			fi
		done

	"""

	stub:
	"""
		touch msa.fa msa_w_controls.fa
	"""
	
}


process GET_COMMON_SNP{
	label "variant_calling"
	publishDir "$outdir/", mode: 'copy'

	input:
	val(variants_timed)  // no touching, helps control timing
	path(path_in)
	val(n_barcodes)
	val(outdir)

	output:
	path("common_snps.vcf.gz"), emit: variants

	script:
	"""
	VARIANTS=\$(find $path_in/ -name "*.gz" -size +0b)
	echo \$VARIANTS

	# Find intersect of all SNPs present in all barcodes
	vcf-isec --apply-filters --force --nfiles =$n_barcodes \$VARIANTS > common_snps.vcf
	bgzip common_snps.vcf

	# Index
	tabix -p vcf common_snps.vcf.gz
	"""

	stub:
	"""
	touch common_snps.vcf.gz
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
			echo ${barcode}.vcf.gz > tmp.txt
			bcftools reheader --samples tmp.txt $path_var/${barcode}.vcf.gz > ${barcode}.vcf.gz
			bcftools index ${barcode}.vcf.gz
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
	cpus=16
	publishDir "$outdir", mode: 'copy'

	input:
	path("msa.fa")
	val(outdir)

	output:
	path("recomb_free.aln"), emit: aln
	path("corrected_tree.nwk"), emit: tree

	script:
	"""
	run_gubbins.py --prefix gubbins --seed 101010 --model GTRGAMMA --recon-model GTRGAMMA \
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
	cpus=16
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
