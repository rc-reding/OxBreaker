// IMPORT MODULES
include { FILTER_READS; FIND_MLST as REFERENCE_MLST; REFERENCE_STATS } from './modules/assembly_utils.nf'
include { EXTRACT_REFERENCE } from './modules/species_screening.nf'

// IMPORT SUB-WORKFLOWS
include { quality_control } from './subworkflows/qc.nf'
include { find_reference } from './subworkflows/species.nf'
include { reference_mapping } from './subworkflows/assembly.nf'
include { variant_caller } from './subworkflows/variant_caller.nf'
include { phylogeny } from './subworkflows/analysis.nf'



// MAIN WORKFLOW
workflow oxbreaker {
	// Channels
		reads = Channel.fromPath("$params.input/*.fastq.gz", checkIfExists: true).map {
				it -> tuple( it.baseName.replace('.fastq',''), it )
				}.toSortedList( a -> a[0] ).flatMap()
		sample = reads.first()
 
	// Main
	main:
		if ( params.model == null & params.clair3 == true ) {
			error "ERROR: You must specify a Clair3 model with '--model'."
		}

		if ( params.min_qual == null ) {
			// No quality filter
			// unless otherwise specified
			min_qual = "0"
		} else {
			min_qual = "$params.min_qual"
		}

		if ( params.min_mq == null ) {
			// Filter reads with mapping quality < 55
			// unless otherwise specified
			min_mq = "55"
		} else {
			min_mq = "$params.min_mq"
		}

		if ( params.min_read_number == null ) {
			// Filter variants with a depth of >=10
			// unless otherwise specified
			min_read_number = "10"
		} else {
			min_read_number = "$params.min_read_number"
		}

		if ( params.min_freq == null ) {
			// Filter variants with a frequency >=90%
			// unless otherwise specified
			min_freq = "0.9"
		} else {
			min_freq = "$params.min_freq"
		}

		// Set reference genome - use one sample only
		if ( !params.reference ) {
			reference = find_reference(sample)

			ref_genbank = reference.genbank
			ref_genome = reference.genome

			REFERENCE_MLST(ref_genome, "$params.output/mlst")

			// Get path only, hacky but works for MLST
			// without breaking anything
			ref_genome = ref_genome.map{ it[1] }
			ref_genbank = ref_genbank.map{ it[1] }
		} else {
			ref_genbank = file("$params.reference")
			reference = EXTRACT_REFERENCE(ref_genbank,
				          	      "$params.output/references")
			ref_genome = reference.genome  //file("$params.output/references/reference.fa")

			REFERENCE_MLST(ref_genome, "$params.output/mlst")
			ref_genome = ref_genome.map{ it[1] }
		}

		// MAIN
		REFERENCE_STATS(ref_genbank,
				"$params.output/references")

		if ( min_qual != "0" ) {
			filtered = FILTER_READS(reads,
						min_qual,
						"$params.output/reads_merged/filtered")
			mapped = reference_mapping(filtered,
						   min_qual,
						   min_mq,
						   ref_genome)
		} else {
			mapped = reference_mapping(reads,
						   min_qual,
						   min_mq,
						   ref_genome)
		}

//		qc = quality_control(reads,
//				     mapped.coverage,
//				     min_qual)

		variant_caller(mapped.asmbl,
			       ref_genome,
			       ref_genbank,
			       mapped.bed,
			       mapped.depth,
			       mapped.depth_asmbl,
			       min_freq,
			       min_read_number,
			       min_mq)
//			       qc.report)		

		phylogeny(variant_caller.out.consensus,
			  variant_caller.out.variants,
			  mapped.coverage,
			  ref_genome)
}
