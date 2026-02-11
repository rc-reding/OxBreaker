// IMPORT MODULES
include { SPECIES_PROFILING_ONT; DOWNLOAD_FASTA } from '../modules/species_screening.nf'



// MAIN WORKFLOW
workflow find_reference {
	take:
		reads
	main:
		SPECIES_PROFILING_ONT(reads,
				      "$params.output/sp_report",
				      "$params.output/sp_report/results")

		reference = DOWNLOAD_FASTA(SPECIES_PROFILING_ONT.out.report,
			     		   "$params.output/references")

	emit:
		genome = reference.genome
		genbank = reference.genbank
}
