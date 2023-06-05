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
	--folder /Users/raul/Downloads/curio_seeker_pipeline_v1.0.3/Mouse_liver/ \
	--sample_name Mouse_liver_R1.fastq.gz \
	--cell_barcodes /Users/raul/Downloads/curio_seeker_pipeline_v1.0.3/Mouse_liver/A0018_040_BeadBarcodes.txt 