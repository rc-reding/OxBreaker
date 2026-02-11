#!/bin/bash

export QT_QPA_PLATFORM=offscreen

DIR="/main/directory/with/data"
echo "Processing ${DIR}..."

# Run oxbreaker
echo "Clearing up cache to save space..."
rm -Rf work

nextflow run main.nf  \
    -profile hpc -resume \
    -entry outbreaker --debug \
    --input $DIR/species/reads \
    # Same numbers if they are not given
    #--min_mq 55 \  
    #--min_freq 0.9 \
    #--min_read_number 10 \
    --phylogeny \
    --output /main/directory/save_here/ \
    # If --reference not given, where is kraken2 DB?
    --db /home/carlos/projects/dbs/kraken/microbial/ \
    # Use this reference (ignores --db)
    --reference /home/carlos/projects/references/s_pasteurianus_WUSP067.gb

echo "Done."
