// IMPORT MODULES
include { MERGE_BARCODES } from './modules/barcode_processing.nf'



// MAIN WORKFLOW
workflow merge_barcodes {
	// Channels
		barcodes = Channel.from("$params.input")
 
	// Main
	main:
		MERGE_BARCODES(barcodes)

}
