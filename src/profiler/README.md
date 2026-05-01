## How the profiler works ##
Credits : Jason Myers
found at sh.d/test-small.sh.d/ART_profiler/illumina in https://github.com/YU-Zhejian/art_modern
___

Based on ARTs simulation model, needs empirical quality score profile derived from real sequencing data. This sequence 
of scripts generates that profile by analysing the distribution of Phred quality scores across read positions from real 
FASTQ files. Scores in Illumina reads vary by position (quality tending to drop towards the end of reads), by nucleotide
base, and by the run. The profiler is intended to capture this variation.

The pipeline consists of 4 scripts that are run through the profiler_illumina.py: 
1. fastqReadAvg.py
    - produces the per-file frequency matrices
2. summation.py
    - merges all the matrices into one combined freq file
3. combinedAvg.py
    - appends weighted average quality scores to the freq file
4. empDist.py
    - converts the freq matrix into a cumulative empirical distribution

Producing a .txt profile file which will be used by the simulator.
Original script allows for usage of multiple threads; mine does not. If that becomes a speed issue with larger datasets
I will need to implement that

### EXAMPLES: ###
	1) create hiseq2k profiles from all *.fq.gz in the directory fastq_dat_dir
		python profiler.py ecoli_hiseq D:\Code\DNA-Sequence-Simulator\src\profiler\fastq_dir fastq
### NOTES: ###
	For paired-end fastq files, e.g., *.fq or *.fq.gz, the filenames of the 1st reads must be *_1.fq/*_1.fq.gz,
	or *.1.fq/*.1.fq.gz, and those of the 2nd reads must be *_2.fq/*_2.fq.gz, or *.2.fq/*.2.fq.gz