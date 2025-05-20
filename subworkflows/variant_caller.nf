include { VARIANT_CALL_CLAIR3; VARIANT_CALL_LONGSHOT } from '../modules/assembly.nf'
include { GENERATE_CONSENSUS } from '../modules/phylogeny.nf'



// MAIN WORKFLOW
workflow variant_caller {
	// Input
	take:
		asmbl
		ref_genome
		depth
		depth_asmbl
		AFthres

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
		variants = VARIANT_CALL_CLAIR3(asmbl.join(depth_asmbl),
					       AFthres, 
					       ref_genome,
					       "$params.output/vcf")
		//variants = VARIANT_CALL_LONGSHOT(asmbl.join(depth_asmbl),
		//				 AFthres,
		//				 ref_genome,
		//				 "$params.output/vcf")

		consensus = GENERATE_CONSENSUS(variants.vcf.join(depth).join(depth_asmbl),
					       asmbl.map{it -> it[1]}, ref_genome,
					       "$params.output/assembly/consensus")

	// Output
	emit:
		variants = variants.vcf
		consensus = consensus.fa
}
