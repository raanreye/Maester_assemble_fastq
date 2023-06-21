'''

This script:
1) Finds the R1 fastq files in the specified folder.
2) Filters for cell barcodes in the cell_barcodes. ...
 ... The cell_barcodes is a file that contains a list of cell barcodes that should be kept.
3) Fastq files are read in. For each read, the code matches the cell barcode to the cell_barcodes.
... If the cell barcode is in the cell_barcodes, the code increments a counter for the corresponding R2 file.
7) A final file.fastq.gz is created with sequences that matched and sequence name is updated to include CellID and UMI.
6) A report containg the number of reads that were filtered by cell_barcodes for each file is outputed.

Input:

  folder = Folder containing fastq files
  sample_name= The names of the files (Ex. file1.fastq.gz file2.fastq.gz), use 'all' if you want to run all the fastqs in your folder.
  cell_barcodes = File (.tsv or .txt) containing cell barcodes

Output:

 file.fastq.gz -> fastq with only sequences from R2 that matched the cell_barcodes. ...
  This final file.fastq.gz will have the barcode (Cell ID or Bead Barcode) ...
  and UMI from R1 and added to the the description (name) of the sequence. 

  file.stats.txt -> summarry of sequences removed based on cell_barcodes.


'''

import subprocess
import os
import sys
import gzip
from Bio import SeqIO, bgzf

def cutf(x, f=1, d="/", **kwargs):
  return sapply(strsplit(x, d), lambda i: i[f], **kwargs)

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def maester_assemble(folder,sample_name,cell_barcodes):

  index_1_2_length = 8

  # Length for CURIO
  #   BB1    linker   BB2   UMI
  # _______ ________ ______ ___ final sequence

  BB1_length = 8
  linker_length = 18
  BB2_length = 6
  curio_umi_length = 9
  polyA = 10

  # Length for 10x genomics
  #   CB       UMI   
  # _______ ________  final sequence

  CB_length = 16
  _10x_umi_length = 12

  # Load file in 1E7 read increments
  batch_size = 1e7


# ____________________________________________________________
# ____________________________________________________________


  # Create folder if it doesn't exist
  if not os.path.exists(folder+'/maester_assemble'):
    os.makedirs(folder+'/maester_assemble')


  # Find R1 fastq files
  r1_files = []
  for root, _, files in os.walk(folder):
    for file in files:

      # Go through all R1.fastq.gz in the folder
      if sample_name == 'all':
        if file.endswith("_R1.fastq.gz"):
          r1_files.append(os.path.join(root, file))

      # Go through only files in sample_name
      elif file in sample_name:
        if file.endswith("_R1.fastq.gz"):
          r1_files.append(os.path.join(root, file))

  # Check if any files.fastq.gz were found
  if r1_files == []:
    raise ValueError("No file.fastq.gz were found.") 

  # Get cell barcodes in the allowlist
  cells = []
  with open(cell_barcodes, "r") as f:
    for line in f:
      cells.append(line.split()[0])

  # Check if any cells were found
  if cells == []:
    raise ValueError("No cells were found.")

  # Process fastq files
  report = {}
  for r1_file in r1_files:
    r2_file = r1_file.replace("_R1", "_R2")
    print("Processing file: {} \n".format(r1_file))

    # Create a list of new records.
    new_records = []

    # Open file.fastq.gz 
    with gzip.open(r1_file, "rt") as f1 , \
         gzip.open(r2_file, "rt") as f2 :


      # Extract cell barcode or bead ID and umi from Read1
      r1_cell = []
      r1_umi = []

      # Subset fastq for 10^7 reads so its faster to run.

      for _, batch in enumerate(batch_iterator(SeqIO.parse(f1, "fastq"), batch_size)):

        # Go through the reads in each batch
        for r1_reads in batch:

          # Check if these reads are from 10x Genomics
          if len(r1_reads) == CB_length + _10x_umi_length:

            # Cell ID
            r1_cell.append(str(r1_reads.seq)[0:CB_length])

            # UMI 
            r1_umi.append(str(r1_reads.seq)[CB_length : CB_length + _10x_umi_length ])


          # Check if these reads are from CURIO
          elif len(r1_reads) == BB1_length+linker_length+BB2_length + curio_umi_length + polyA:

            # Bead ID 1 & 2
            BB_1 = str(r1_reads.seq)[0:BB1_length]
            BB_2 = str(r1_reads.seq)[BB1_length+linker_length:BB1_length+linker_length+BB2_length]
            r1_cell.append(BB_1 + BB_2)

            # UMI 
            r1_umi.append(str(r1_reads.seq)[BB1_length+linker_length+BB2_length:BB1_length+linker_length+BB2_length + curio_umi_length ])

          else:
            raise ValueError("Length of R1 reads didn't match 10x Genomics or CURIO.")
        
      # What cells match the whitelist
      set_cells = set(cells)
      index_of_cell = set([i for i, item in enumerate(r1_cell) if item in set_cells])

      cnt = 0
      # Subset fastq for 10^7 reads so its faster to run.
      for _, batch in enumerate(batch_iterator(SeqIO.parse(f2, "fastq"), batch_size)):

        # Open batch R2 and subset it based of cells that matched and change description ...
        # to include cell ID and UMI
        for r2_reads in batch:

          # Only keep cells that were matched
          if cnt in index_of_cell:

              hold_description = r2_reads.description.split(' ')[0] + '_' + r2_reads.description.split(' ')[1][(-1-2*index_1_2_length):]+ '_' + r1_cell[cnt] + '_' +r1_umi[cnt]
              r2_reads.id = hold_description
              r2_reads.description = '<unknown description>' # standard for SeqIO
              
              # Add the record to the list of new records.
              new_records.append(r2_reads)
            
          cnt += 1


      # Write the modified.fastq.gz
      with bgzf.BgzfWriter(folder+'/maester_assemble/' + args.sample_name.split('.')[0][0:-3] + '.fastq.gz', "wb") as output_handle:
        SeqIO.write(sequences=new_records, handle=output_handle, format="fastq")

      # Write report
      with open(folder+'/maester_assemble/' + args.sample_name.split('.')[0][0:-3] + ".stats.txt", "w") as f:
        f.write("all\tfiltered\tfraction\n")
        f.write("{}\t{}\t{:.2f}\n".format(len(r1_cell), len(index_of_cell), len(index_of_cell) / len(r1_cell)))

if __name__ == "__main__":

  import argparse

  # Get arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--folder",type=str, help="Folder containing fastq files")
  parser.add_argument("--sample_name",type=str, help="Sample name")
  parser.add_argument("--cell_barcodes",type=str, help="File (.tsv) containing cell barcodes")

  args = parser.parse_args()

  maester_assemble(args.folder,args.sample_name,args.cell_barcodes)







