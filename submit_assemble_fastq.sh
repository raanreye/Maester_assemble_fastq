#!/bin/bash
#submit using: bsub < /home/raanreye/Submit_FastqFragment.sh
#BSUB -J Well_job
#BSUB -o Well_job.%J.out
#BSUB -e Well_job.%J.error
#BSUB -M 200000 # Reqesting 10GB RAM
#BSUB -R "span[hosts=1] rusage [mem=200000]"
#BSUB -q normal #input gpu if using gpu

module load python/3.9.1

python /Users/raul/Documents/GitHub/Maester_assemble_fastq/Maester_assemble_fastq.py \
	--folder /Volumes/GoogleDrive/.shortcut-targets-by-id/1CzHSa6eWu2vP_5qpBUKrCi8hMOj0gF1C/gier_shared/projects/gej_model/raul/subset_test/ \
	--sample_name subset_R1 \
	--cell_barcodes /Volumes/GoogleDrive/.shortcut-targets-by-id/1CzHSa6eWu2vP_5qpBUKrCi8hMOj0gF1C/gier_shared/projects/gej_model/raul/barcodes.tsv \
	--BB1_length 8 \
	--linker_length 18 \
	--BB2_length 6 \
	--umi_length 4