
DIR="/mnt/baradur/carlos/outbreaks/Dorado_v1.0.2/mrsa_R10"
echo "Processing ${DIR}..."
	echo "Processing $DIR"
	nextflow run preprocess.nf $1 \
		-entry merge_barcodes \
		--input $DIR/reads/ \
		--output $DIR/reads_merged
