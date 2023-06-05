#!/bin/bash
#submit using: bsub < /home/raanreye/Submit_FastqFragment.sh
#BSUB -J Well_job
#BSUB -o Well_job.%J.out
#BSUB -e Well_job.%J.error
#BSUB -M 200000 # Reqesting 10GB RAM
#BSUB -R "span[hosts=1] rusage [mem=200000]"
#BSUB -q normal #input gpu if using gpu

module load python/3.9.1

python /path/to/Maester_assemble_fastq.py \
	--folder /path/to/fastq/folder/ \
	--sample_name file1_R1.fastq.gz file2_R1.fastq.gz \
	--cell_barcodes /path/to/BeadBarcodes.txt 