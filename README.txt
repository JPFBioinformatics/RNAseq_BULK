Pipeline to process bulk RNA sequencing data
built to be ran on HPC clusers with SLURM job handling for efficient 

Works for 2 read illumina sequencing
built for use on HPC clusters or on local machines, but uses STAR which is memory-intensive so without a lot of spare memory you might not want to run locally

Just run run_pipeline.sh to enter in at main.py and run sequencing as it is automatically setup

=========================================================== Naming Conventions ===========================================================
    DO NOT:
        include any "." in your filename except for on file extenstions
            sample1.v1.fastq.gz is wrong but sample1_v1.fastq.gz works well
    DO:
        make sure to differentiate gropued fastq files using _R1 and _R2
            sample1_v1_R1.fastq.gz and sample1_v1_R2.fastq.gz will be grouped together
            this allows for automatic sample gropu detection and facilitates efficeint SLURM usage
    
    We do allow for the user to specify two files to group together, but it is clunker for SLURM automation so it is meant to be used sparingly

============================================================= main arguments =============================================================

--root                  Used to specify the location of the root folder (RNAseq_BULK) to facilitate good HPC cluster use
--indir                 Path object that points to the directory where your .fastq or .fastq.gz (preferred) files are stored
--sample1               filename for sample 1 forward read, optional.  If not specified then the script will process all .fastq/.fastq.gz files
--sample2               located in at --indir.  Specify only --sample1 to take advantage of automatic forward/reverse reading for SLURM HPC automation
--steps                 "trim", "align", or "count" specifies the part of the pipeline you want to run for easy HPC use if not specified then all steps will run
--cleanup               you REALLY want to include this in your CLI call.  It makes it so that intermediate files are automatically deleted after used,
                        saving a lot of memory space.  If you do --cleanup False then it will save the intermediate files (BAM, BAI etc...) to the root
                        folder so you can access them later
