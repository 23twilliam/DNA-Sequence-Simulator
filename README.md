# DNA-Sequence-Simulator
Major Project wim35

The simulator is run using python(3) main.pu run

With optional flags including  
-i, --input, for the input genome file path  
-p, --plaform, for the desired platform to be simulated  
-l, --read_length, for the desired length of the reads  
-c, --coverage, to set the desired coverage depth  
-f, --format, to select between output formats (FASTA, FASTQ)  
-o, --output, for the desired output file path  
--profile, if you would like to use a custom error profile

In built defaults make all of these optional.

Example commands:  
python main.py run -i genome.fasta  

python main.py run -c 20 -l 150 -p pacbio  

python main.py -f fasta -o some_dir -p nanopore
