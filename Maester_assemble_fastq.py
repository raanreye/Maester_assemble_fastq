'''

This script:
1) Finds the R1 fastq files in the specified folder.
2) Filters for cell barcodes in the allowlist. The allowlist is a file that contains a list of cell barcodes that should be kept.
3) Fastq files are reads in 10 million read increments. 
... For each read, the code matches the cell barcode to the allowlist.
... If the cell barcode is in the allowlist, the code increments a counter for the corresponding R1 file.
6) A report containg the number of reads that were filtered by cell barcode for each R1 file is outputed.



Input:



Output:



'''

import subprocess
import os
import sys
import Bio
import gzip
from Bio import SeqIO

def cutf(x, f=1, d="/", **kwargs):
  return sapply(strsplit(x, d), lambda i: i[f], **kwargs)

def maester_assemble(folder,sample_name,cell_barcodes):

  # Length for CURIO
  #   BB1    linker   BB2   UMI
  # _______ ________ ______ ___ final sequence

  BB1_length = 8
  linker_length = 18
  BB2_length = 6
  curio_umi_length = 4

  # Length for 10x genomics
  #   CB       UMI   
  # _______ ________  final sequence

  CB_length = 16
  _10x_umi_length = 12





  # Find R1 fastq files
  r1_files = []
  for root, _, files in os.walk(folder):
    for file in files:
      if file.endswith("_R1.fastq.gz"):
        r1_files.append(os.path.join(root, file))

  # Get cell barcodes in the allowlist
  cells = []
  with open(cell_barcodes, "r") as f:
    for line in f:
      cells.append(line.strip())

  # Check if any cells were found
  if cells == []:
    raise ValueError("No cells were found.")

  # Process fastq files
  report = {}
  for r1_file in r1_files:
    r2_file = r1_file.replace("_R1", "_R2")
    print("Processing file: {}".format(r1_file))

    # Load file in 1E7 read increments
    with gzip.open(r1_file, "rt") as f1 , \
         gzip.open(r2_file, "rt") as f2 :

      # Extract cell barcode or bead ID and umi from Read1
      r1_cell = []
      r1_umi = []
      for r1_reads in SeqIO.parse(f1, "fastq"):

        # Check if these reads are from 10x Genomics
        if len(r1_reads) == CB_length + _10x_umi_length:

          # Cell ID
          r1_cell.append(str(r1_reads.seq)[0:CB_length])

          # UMI 
          r1_umi.append(str(r1_reads.seq)[CB_length : CB_length + _10x_umi_length ])


        # Check if these reads are from 10x Genomics
        elif len(r1_reads) == BB1_length+linker_length+BB2_length + curio_umi_length:
          # Bead ID 1 & 2
          BB_1 = str(r1_reads.seq)[0:BB1_length]
          BB_2 = str(r1_reads.seq)[BB1_length+linker_length:BB1_length+linker_length+BB2_length]
          r1_cell.append(BB_1 + BB_2)

          # UMI 
          r1_umi.append(str(r1_reads.seq)[BB1_length+linker_length+BB2_length:BB1_length+linker_length+BB2_length + curio_umi_length ])

        else:
          raise ValueError("Length of R1 reads didn't match 10x Genomics or CURIO.")




      r2_array = []
      r2_umi = []
      for r2_reads in SeqIO.parse(f2, "fastq"):
        r2_array.append(str(r2_reads.seq))
        r2_umi.append(str(r2_reads.seq))


      cells_matched = list(set(r1_cell) - set(cells))

      # print(cells_matched)

      # output_handle = open("DIFF.fastq","w")

      # SeqIO.write(needed_reads,output_handle,"fastq")

      # output_handle.close()

  # # Write report
  # with open(args.sample_name + ".stats.txt", "w") as f:
  #   f.write("all\tfiltered\tfraction\n")
  #   for r1_file, count in report.items():
  #     f.write("{}\t{}\t{:.2f}\n".format(count, count, count / 10**7))

if __name__ == "__main__":

  import argparse

  # Get arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--folder",type=str, help="Folder containing fastq files")
  parser.add_argument("--sample_name",type=str, help="Sample name")
  parser.add_argument("--cell_barcodes",type=str, help="File (.tsv) containing cell barcodes")

  args = parser.parse_args()

  maester_assemble(args.folder,args.sample_name,args.cell_barcodes)







