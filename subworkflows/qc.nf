include { QC_ONT } from '../modules/assembly_utils.nf'
include { PLOT_QC; PLOT_QC as PLOT_FILTERED_QC } from '../modules/plot_figures.nf'
// Cannot call a process more than once, rename to re-use
include {QC_ONT as FILTERED_QC_ONT} from '../modules/assembly_utils.nf'



// MAIN WORKFLOW
workflow quality_control{
	// Function
	take:
		reads
		coverage
		min_qual

	main:
		control = QC_ONT(reads,
				 "$params.output/qc")

		PLOT_QC("$params.output/qc",
			control.qc.collect(),
			coverage.collect(),
			"$params.output/qc")

	emit:
		report = control.qc
}
