// IMPORT MODULES
include { FIND_MLST as REFERENCE_MLST } from './modules/assembly_utils.nf'

// IMPORT SUB-WORKFLOWS
include { quality_control } from './subworkflows/qc.nf'
include { find_reference } from './subworkflows/species.nf'
include { reference_mapping; denovo_assembly } from './subworkflows/assembly.nf'
include { variant_caller } from './subworkflows/variant_caller.nf'
include { phylogeny } from './subworkflows/analysis.nf'



// MAIN WORKFLOW
workflow outbreaker {
	// Channels
		reads = Channel.fromPath("$params.input/reads_merged/*.fastq.gz", checkIfExists: true).map {
				it -> tuple( it.baseName.replace('.fastq',''), it )
				}.toSortedList( a -> a[0] ).flatMap()
		sample = Channel.fromPath("$params.input/reads_merged/*.fastq.gz", checkIfExists: true).map {
				it -> tuple( it.baseName.replace('.fastq',''), it )
				}.first()
 
	// Main
	main:

		if ( params.min_qual == null) {
			// No filtering by default
			// unless otherwise specified
			min_qual = "0"
		} else {
			min_qual = "$params.min_qual"
		}

		/* INTERESTING: positives seem to decline with 
		 * higher min_freq up to ~50-60%. However,
		 * higher min_freq requirement results
		 * in a modest increase in number of SNP.
		 * But it is not consistent---depends on sp.
		 * More robust is a dip occurring at 90% freq.
		 * Information is being _thrown_ at min_freq
		 * above 50% or 90%? 
		 */

		if ( params.min_freq == null ) {
			min_freq = "0.95"
		} else {
			min_freq = "$params.min_freq"
		}

		// Set reference genome - use one sample only
		if ( !params.reference ) {
			reference = find_reference(sample)
			ref_genome = reference.genome

			REFERENCE_MLST(ref_genome, "$params.output/mlst")

			// Get path only, hacky but works for MLST
			// without breaking anything
			ref_genome = ref_genome.map{ it[1] }
		} else {
			ref_genome = file("$params.reference")
			REFERENCE_MLST(tuple("reference", ref_genome),
				       "$params.output/mlst")
		}

		// MAIN
		if ( "$params.whole_genome" == true ) {
			dn = denovo_assembly(reads)

			// Quality control
			quality_control(reads,
					dn.coverage,
					min_qual)
	
			// TODO: ue minimap2 / something to sort contig based on references
			variants = ''
		} else {
			mapped = reference_mapping(reads,
						   ref_genome)

			quality_control(reads,
					mapped.coverage,
					min_qual)
	
			variant_caller(mapped.asmbl,
				       ref_genome,
				       mapped.depth,
				       mapped.depth_asmbl,
				       min_freq)
		}
		

		phylogeny(variant_caller.out.consensus,
			  variant_caller.out.variants,
			  mapped.coverage,
			  ref_genome)
}
