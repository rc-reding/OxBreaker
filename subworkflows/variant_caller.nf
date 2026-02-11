include { VARIANT_CALL_CLAIR3; VARIANT_CALL_BCFTOOLS } from '../modules/assembly.nf'
include { RETRIEVE_CONTEXT } from '../modules/assembly_utils.nf'
include { GENERATE_CONSENSUS } from '../modules/phylogeny.nf'
include { PLOT_SNP_DEPTH_DISTR_CLAIR3; PLOT_SNP_DEPTH_DISTR_BCFTOOLS; PLOT_SNP_QUAL_DISTR } from '../modules/plot_figures.nf'



// MAIN WORKFLOW
workflow variant_caller {
	// Input
	take:
		asmbl
		ref_genome
		ref_genbank
		bed_file
		depth
		depth_asmbl
		min_freq
		min_read_number
		min_mapping_score
//		qc_status

	// Main
	main:
		/* BEWARE https://github.com/HKU-BAL/Clair3/issues/123#issuecomment-1185498217
                 * Clair3 does NOT filter SNPs based on frequency, SNP length, etc given
                 * by the user.
                 *
                 * Rather, these settings filter _candidates_ that are fed to the deep
                 * learning model. The final calls do not necessarily meet the settings
                 * given by the user.
                 */
	if ( params.clair3 == true ) {
		variants = VARIANT_CALL_CLAIR3(asmbl.join(depth_asmbl),
					       min_freq,
					       min_read_number,
					       min_mapping_score,
					       ref_genome,
					       "$params.output/vcf")

		PLOT_SNP_DEPTH_DISTR_CLAIR3(variants.vcf_report,
				     depth.collect(),
				    "$params.output/depth",
				    "$params.output/vcf/figures")
	} else {
		variants = VARIANT_CALL_BCFTOOLS(asmbl.join(depth_asmbl),
					       min_freq,
					       min_read_number,
					       min_mapping_score,
					       ref_genome,
					       bed_file,
					       "$params.output/vcf")

		PLOT_SNP_DEPTH_DISTR_BCFTOOLS(variants.vcf_report,
				     depth.collect(),
				    "$params.output/depth",
				    "$params.output/vcf/figures")
	}

		// REQUIRES MODIFIED VERSION OF NANOPLOT TO EXPORT STD
//		PLOT_SNP_QUAL_DISTR(variants.vcf,
//				    qc_status.collect(),
//				    "$params.output/qc",
//				    "$params.output/vcf/figures")

		variants4msa = variants.vcf.join(asmbl).groupTuple().map{
							barcode, asmbl, vcfs -> tuple(barcode,
										      asmbl[0],
										      vcfs[0])
							}
		
		consensus = GENERATE_CONSENSUS(variants4msa.join(depth).join(depth_asmbl),
					       ref_genome,
					       bed_file,
					       min_read_number,
					       "$params.output/assembly/consensus")

		variants4ctx = variants.vcf_report.join(asmbl).groupTuple().map{
							barcode, asmbl, vcfs -> tuple(barcode,
										      asmbl[0],
										      vcfs[0])
							}

		compiled_input = variants4ctx.join(consensus.fa).groupTuple().map{
				barcode, vcfs, asmbl, consensus_fa -> tuple(barcode,
									vcfs[0],
									asmbl[0],
									consensus_fa[0])
			}
		
		RETRIEVE_CONTEXT(compiled_input,
				 ref_genbank,
				 "$params.output/vcf")

	// Output
	emit:
		variants = variants.vcf_report
		consensus = consensus.fa
}
