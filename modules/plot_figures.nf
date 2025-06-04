process PLOT_PHYLOGENY {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	val(phylogeny_dir)
	path(phy_tree)  // Only to ensure the process runs after phylogeny
	val(mlst_dir)
	val(outdir)

	output:
	path("phylogeny.pdf"), emit: phylogeny

	script:
	"""
	python3 $params.bin/plot_dendrogram.py $phylogeny_dir/ ./ $mlst_dir/
	"""

	stub:
	"""
	touch phylogeny.pdf
	"""
}


process PLOT_SNP_DISTANCES {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	path(snp_dists) 
	val(outdir)

	output:
	path("snp_histogram.pdf"), emit: snp_hist
	path("distance_matrix.pdf"), emit: dist_matrix

	script:
	"""
	python3 $params.bin/plot_snp_histogram.py $snp_dists $outdir
	"""

	stub:
	"""
	touch snp_histogram.pdf distance_matrix.pdf
	"""
}


process PLOT_SNP_DEPTH_DIST {
	label "phylogeny_gubbins"
	publishDir "$outdir", mode: 'copy'

	input:
	tuple val(barcode), path(vcf_file)
	path(depth_dir)
	val(outdir)

	output:
	path("${barcode}_snp_dist.pdf"), emit: depth_dist

	script:
	"""
	python3 $params.bin/plot_snp_depth_dist.py $vcf_file $depth_dir/${barcode}_depth.csv
	"""

	stub:
	"""
	touch ${barcode}_snp_dist.pdf
	"""
}



process PLOT_COVERAGE {
	label "phylogeny"
	publishDir "$outdir", mode: 'copy'

	input:
	val(depth_dir)
	val(coverage)  // Only to ensure process runs after depth
	val(outdir)

	output:
	path("samples_breadth_coverage.pdf"), emit: breadth
	path("samples_coverage_depth.pdf"), emit: depth

	script:
	"""
	python3 $params.bin/plot_coverage.py $depth_dir/ ./
	"""

	stub:
	"""
	touch samples_coverage_depth.pdf samples_breadth_coverage.pdf
	"""
}

process PLOT_QC {
	label "phylogeny"
	publishDir "$outdir", mode: 'copy'

	input:
	val(qc_dir)
	val(qc_data)  // Only to ensure process runs after qc complete
	val(coverage)  // Only to ensure process runs after depth complete
	val(outdir)

	output:
	path("samples_qc.pdf"), emit: breadth

	script:
	"""
	python3 $params.bin/plot_qc.py $qc_dir/ ./
	"""

	stub:
	"""
	touch samples_qc.pdf
	"""
}
