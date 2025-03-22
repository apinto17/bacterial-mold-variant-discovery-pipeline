## Overview

This is a pipeline that is used to identify and analyze novel sequences of mold strains. The strains are then aligned to the reference sequence for variant discovery. 

- dgorgon_reference.fa - The wildtype reference sequence of D.gorgon
- harrington_clinical_data.txt - The clinical data of the samples including a unique patient name, the color of the mold, and a sequence barcode.
- hawkins_pooled_sequences.fastq - The fastq data that contains the sequences of all the samples pooled together.


## Pipeline steps (all in pipeline.py)

1. Reads in the pooled fastq and outputs 50 new fastq files (for 50 patients) that contain only sequences belonging to that sample. The samples are identified using the barcode and trimmed of poor reads according to the fastq file.

2. Perform alignment on each FASTQ to the reference sequence. The align function calls the bwa command. This will transform the FASTQs to samfiles.

3. Compress alignments. Convert the samfiles to bamfiles. Uses samtools view to convert the SAM file into a BAM file. BAM is the binary format corresponding to the SAM text format.

4. Pileup sequences and identify SNP's. Uses the python pysam pileup function to discover variants in each sorted bam file. 

5. Create a report that outputs what nucleotide position and mutation is responsible for each color of the mold. Also prints out the number of sequences that were used for each sample. 

ex. Sample Tim had a green mold, 320 reads, and had 32% of the reads at position 23 had the mutation C.

## Dependencies

- make sure you have the hawkins_pooled_sequences.fastq file, or the file you're interested in, at the same directory level as pipeline.py

- make sure you have "helper_scripts" at the same directory level

- make sure you have pysam installed

## to Run

`python pipeline.py -f hawkins_pooled_sequences.fastq`