// IMPORT MODULES
include { CONSENSUS_2_MSA; GET_COMMON_SNP; PREALIGN_GENOMES; GET_SAMPLES_NUMBER;
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
		if ( "$params.whole_genome" == true ) {
			dist = DISTANCE_MATRIX(wgMLST_msa.collect(), 
					       "$params.output/phylogeny/wg")
							
			dist_wo_controls = DISTANCE_MATRIX_wo_CTRL(wgMLST_msa_wo_controls.collect(),
								   "$params.output/phylogeny/wg/no_controls")
							
			// Phylogeny
			phy = GENERATE_PHYLOGENY_GUBBINS(wgMLST_msa.collect(),
							 "$params.output/phylogeny/wg")
			// Maximum-Likelihood (ML) method with RaxML
			//phy_ML = GENERATE_PHYLOGENY(wgMLST_msa.collect(), "$params.output/phylogeny/wg")
			//phy = CORRECT_RECOMBINATION(phy_ML.tree, wgMLST_msa.collect(),
			//						    "$params.output/phylogeny/wg")
			
			// SNP distance required for plotting, wait for completion
			dist.snp.collect()
			dist_wo_controls.snp.collect()

			// Plot
			PLOT_PHYLOGENY("$params.output/phylogeny/wg",
				       phy.tree,
				       "$params.output/mlst/wg",
				       "$params.output/phylogeny/wg")
		} else {
		  	rehead_variants = PREALIGN_GENOMES(coverage,
							   variants.collect(),
							   "$params.output/vcf",
							   "$params.output/vcf/rehead")
				
			// Alignment must wait for all reads to be pre-processed
			barcodes = rehead_variants.preprocessed.map{ it -> it[1] }.collect()
			n_barcodes = barcodes.size()

			// Get variants common to all samples
			common =  GET_COMMON_SNP(barcodes,
						 "$params.output/vcf/rehead",
						 n_barcodes,
						 "$params.output/vcf")
		
			msa = CONSENSUS_2_MSA(consensus.collect(),
					      "$params.output/assembly/consensus",
					      ref_genome,
					      "$params.output/alignments")

			dist = DISTANCE_MATRIX(msa.alignment.collect(),
					       "$params.output/phylogeny")

			// Phylogeny
//			phy = GENERATE_PHYLOGENY(msa.alignment.collect(),
//						 "$params.output/phylogeny")

			// Maximum-Likelihood (ML) method with RaxML
			//phy_ML = GENERATE_PHYLOGENY_RAxML(msa.alignment.collect(),
			//			    "$params.output/phylogeny")
			//
			//phy = CORRECT_RECOMBINATION(phy_ML.tree,
			//			    msa.alignment.collect(),
			//			    "$params.output/phylogeny")
			
			// SNP distance required for plotting, wait for completion
//			dist.snp.collect()
			
			// Plots
			PLOT_SNP_DISTANCES(dist.snp,
					   "$params.output")

//			PLOT_PHYLOGENY("$params.output/phylogeny",
//				       phy.tree,
//				       "$params.output/mlst",
//				       "$params.output")

			if ( params.debug == true ) {
				// Re-do phylogeny, SNP distance comparison, and figures
				// including reference genome and controls where available
				// for sanity checks
				dbug_dist = DEBUG_DISTANCE_MATRIX(msa.alignment_w_controls.collect(),
								  "$params.output/phylogeny/debug")
							
				PLOT_DISTANCES_DEBUG(dbug_dist.snp,
						     "$params.output/phylogeny/debug")

				dbug_phy = DEBUG_PHYLOGENY(msa.alignment_w_controls.collect(),
							   "$params.output/phylogeny/debug")
				PLOT_PHYLOGENY_DEBUG("$params.output/phylogeny/debug",
						     dbug_phy.tree,
						     "$params.output/mlst",
						     "$params.output/phylogeny/debug")
			}
		}
}
