#!/bin/bash
conda activate telescope
python cellsplit_sc.py --bam /efs/projects/align-pbmc/results/pbmc_10k_v3_fastqs/pbmc_10k_v3_GDC38.Aligned.sortedByCoord.out.bam --brcds /efs/projects/align-pbmc/results/pbmc_10k_v3_fastqs/pbmc_10k_v3_GDC38.Solo.out/Gene/filtered/barcodes.tsv --dest /fsx/users/santiago2/cellsplit_sc_results/ --log /fsx/users/santiago2/cellsplit_log.tsv
