include { MAP_REFERENCE } from '../modules/assembly.nf'
include { CALCULATE_COVERAGE_DEPTH; BED_FROM_BLAST } from '../modules/assembly_utils.nf'
include { PLOT_COVERAGE; PLOT_DEPTH_DISTR } from '../modules/plot_figures.nf'



workflow reference_mapping {
	take:
		reads
		min_qual
		min_mapping_score
		ref_genome

	main:
		bed = BED_FROM_BLAST(ref_genome,
				   "$params.output/assembly")


		mapped_asmbl = MAP_REFERENCE(reads,
					     ref_genome,
					     bed.file,
					     min_mapping_score,
					     "$params.output/assembly")

		barcode_depth = CALCULATE_COVERAGE_DEPTH(mapped_asmbl.bam,
							 min_mapping_score,
							 "$params.output/depth")

		PLOT_COVERAGE("$params.output/depth",
			      barcode_depth.cov.collect(),
			      "$params.output/depth")

		PLOT_DEPTH_DISTR(barcode_depth.tsv,
			         "$params.output/depth/figures")

	emit:
		asmbl = mapped_asmbl.bam
		coverage = barcode_depth.cov
		depth = barcode_depth.tsv
		depth_asmbl = barcode_depth.depth
		bed = bed.file
}

