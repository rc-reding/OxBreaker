// IMPORT MODULES
include { CONSENSUS_2_MSA; PREALIGN_GENOMES; GET_SAMPLES_NUMBER;
	  DISTANCE_MATRIX; DISTANCE_MATRIX as DEBUG_DISTANCE_MATRIX;
	  GENERATE_PHYLOGENY; GENERATE_PHYLOGENY as DEBUG_PHYLOGENY } from '../modules/phylogeny.nf'
include { PLOT_PHYLOGENY; PLOT_PHYLOGENY as PLOT_PHYLOGENY_DEBUG; 
	  PLOT_SNP_DISTANCES; PLOT_SNP_DISTANCES as PLOT_DISTANCES_DEBUG } from '../modules/plot_figures.nf'

// IMPORT PROCESSES
include { MLST } from './mlst.nf'



// MAIN WORKFLOW
workflow phylogeny {
	take:
		consensus
		variants
		coverage
		ref_genome

	main:
		MLST = MLST(consensus, ref_genome)

		// Gubbins allows SH-test, ClonalFrameML does not.
		// Length correction similar otherwise, with Gubbins being
		// a bit more conservative than ClonalFrameML
		rehead_variants = PREALIGN_GENOMES(coverage,
						   variants.collect(),
						   "$params.output/vcf",
						   "$params.output/vcf/rehead")
			
		// Alignment must wait for all reads to be pre-processed
		barcodes = rehead_variants.preprocessed.map{ it -> it[1] }.collect()
		n_barcodes = barcodes.size()

		msa = CONSENSUS_2_MSA(consensus.collect(),
				      "$params.output/assembly/consensus",
				      ref_genome,
				      "$params.output/alignments")

		dist = DISTANCE_MATRIX(msa.alignment.collect(),
				       "$params.output/phylogeny")

		if ( params.phylogeny == true ) {
			// Phylogeny
			phy = GENERATE_PHYLOGENY(msa.alignment.collect(),
						 "$params.output/phylogeny")

			PLOT_PHYLOGENY("$params.output/phylogeny",
				       phy.tree,
				       "$params.output/mlst",
				       "$params.output")
		}

		// SNP distance required for plotting, wait for completion
		dist.snp.collect()
		
		// Plots
		PLOT_SNP_DISTANCES(dist.snp,
				   "$params.output")


		if ( params.debug == true ) {
			// Re-do phylogeny, SNP distance comparison, and figures
			// including reference genome and controls where available
			// for sanity checks
			dbug_dist = DEBUG_DISTANCE_MATRIX(msa.alignment_w_controls.collect(),
							  "$params.output/phylogeny/debug")			

			// SNP distance required for plotting, wait for completion
			dbug_dist.snp.collect()

			PLOT_DISTANCES_DEBUG(dbug_dist.snp,
					     "$params.output/phylogeny/debug")
		
			if ( params.phylogeny == true ) {
				dbug_phy = DEBUG_PHYLOGENY(msa.alignment_w_controls.collect(),
							   "$params.output/phylogeny/debug")

				PLOT_PHYLOGENY_DEBUG("$params.output/phylogeny/debug",
						     dbug_phy.tree,
						     "$params.output/mlst",
						     "$params.output/phylogeny/debug")
			}
		}
}
