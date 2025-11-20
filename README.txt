====================================================== Index =====================================================

Sections:
    General                         General use information, basics of what the pipeline does and how to run with defaults
    Naming Conventions              what to do/not do when naming your input fastq files (very important for automation)
    Arguments                       description of the arguments that you need to pass with the script call(s) in your command line
    Config                          information on config.yaml, how you can get the pipeline to run differently if you don't like the defaults, probably don't touch
                                    config.yaml if you don't know a lot about the tools used in the pipeline
    Testing                         where the test files in this repo came from and how to set up testing for the pipeline

===================================================== General ====================================================

Pipeline to process bulk RNA sequencing data
Built to be ran on HPC clusers with SLURM job handling for efficient processing, but can be ran locally as well, just keep in mind that STAR aligner uses a lot
of memory so you need lots of extra memory and this will take a long time if you don't have a lot of processing power
Built for 2 read illumina sequencing
PLEASE consult Naming Conventions section of readme at least before running

You must also either place your genome fasta file and GTF file in the reference dir  of this pipeline (RNAseq_BULK), withe their names being set to "genome.fa" and "genes.gtf" 
respectively, or specify the dir that these files are in and their names in the config.yaml before running the pipeline

to run this pipeline with default parameters here are some example commands for different use cases

        NOTES:

        root folder is the root of this pipeline, RNAseq_BULK or whatever you renamed it
        ALL paths must be ABSOLUTE (from filesystem root)

    1. running the pipeline on a set of samples on an HPC cluster

python launcher.py \
    --root /aboslute/path/to/root/folder \
    --indir /home/yourname/RNAseq_bulk/raw_data \
    --steps trim align count

    2. running the pipeline locally on a set of samples (not recommended unless you have a lot of memory/processing power)

python main.py \
    --root /Users/john/RNAseq_bulk \
    --indir /Users/john/RNAseq_bulk/raw_data \
    --steps trim align count

    3. running the pipeline on a single set of samples you specify
            if you follow naming conventions then you can just specify --sample1 and omit --sample2 line

python main.py \
    --root /Users/john/RNAseq_bulk \
    --indir /Users/john/RNAseq_bulk/raw_data \
    --sample1 sampleA_S1_L001_R1_001.fastq.gz \
    --sample2 sampleA_S1_L001_R2_001.fastq.gz \
    --steps trim align count

for any of these 

================================================ Naming Conventions ================================================

    DO NOT:
        include any "." in your filename except for on file extenstions
            sample1.v1.fastq.gz is wrong but sample1_v1.fastq.gz works well
    DO:
        make sure to differentiate gropued fastq files using _R1 and _R2
            sample1_v1_R1.fastq.gz and sample1_v1_R2.fastq.gz will be grouped together
            this allows for automatic sample gropu detection and facilitates efficeint SLURM usage

    We do allow for the user to specify two files to group together, but it is clunker for SLURM automation so it is meant to be used sparingly

===================================================== Arguments ===================================================

launcher.py

    python script that will configure the SLURM args based on the input data/config.yaml file parameters

        --root                  root arg to pass to main.py
        --indir                 indir arg to pass to main.py
        --steps                 steps arg to pass to main.py
        --runScript             slurm (.sh) script that launcher is building

run_pipeline.sh

    SLURM script that handles resource allocation for HPC cluster, directing each available node to process a single sample set, works for local single node machines as well
    no args itself, launcher.py handles automatic configuration, but this script will detect/create a scratch/temp dir for intermediate files

main.py

    This script runs the actual pipeline, uses FastP for QC/trimming, STAR for index generation/aligning, samtools to sort/filter/index aligned bam files, and featureCounts
    to generate count text files for downstream analysis

        --root                  Used to specify the location of the root folder (RNAseq_BULK) to facilitate good HPC cluster use
        --indir                 Path object that points to the directory where your .fastq or .fastq.gz (preferred) files are stored
        --sample1               filename for sample 1 forward read, optional.  If not specified then the script will process all .fastq/.fastq.gz files
        --sample2               located in at --indir.  Specify only --sample1 to take advantage of automatic forward/reverse reading for SLURM HPC automation
        --steps                 "trim", "align", and/or "count" specifies the part of the pipeline you want to run for easy HPC use if not specified then all steps will run

=================================================== Config file ====================================================

This pipeline is built to work off of the config.yaml file located in the root folder
you can edit this file as directed within the file, but here is a more in-depth description of what each field in config means
keep in mind for booleans yaml is canse insensitive for ONLY the FIRST CHARACTER (so true = True and false = False), but other than that will not recognise as booleans 
    tRue =/ True and faLsE =/ False for Example
    we have a check_bool called in launcher.py for this, but not in main so if you manually call main be aware of this, or add a cfg.check_bool() call to the beginning of
    main if you are manually running this
Default values are already set in the config so don't mess witht them unless you know what you are doing

Project Configuration (-project:)

    -name
    Name of the RNA-seq run. A directory with this name will be created under runs/, containing data/, logs/, and reference/ as well as ouptus, must be unique for each run.

    -save_files
    Whether to keep intermediate BAM/CRAM files after producing counts.
        True = keep BAM/CRAM files
        False = delete them after featureCounts completes

    -save_type
    Which file format to save when save_files is True.
        "bam" = save BAM files
        "cram" = save compressed CRAM files (requires genome FASTA for decoding)

Reference Configuration (-reference:)

    -ref_dir
    Directory under the project folder where reference files are stored.

    -genome_fasta
    Reference genome FASTA used for:
        1. STAR index building
        2. STAR alignment
        3. CRAM creation/decoding
        4. featureCounts
        Must be uncompressed (.fa or .fasta).

    -gtf_file
    Genome annotation GTF file used by featureCounts to assign reads to genes.

    -star_index
    Directory containing STAR’s genome index. If empty, pipeline builds it automatically.

Tool Configuration (-tools:)

    FastP Configuration (-tools.fastp)

        -length_required
        Minimum read length after trimming. Reads shorter than this are removed.

        -qualified_quality_phred
        Phred score threshold for trimming low-quality bases from read ends.

        -specify_adapter
        False = fastp automatically detects Illumina adapters.
        True = manually specify adapter_sequence and adapter_sequence_r2.

        -adapter_sequence
        Manual adapter sequence for R1 reads (only used if specify_adapter = True).

        -adapter_sequence_r2
        Manual adapter sequence for R2 reads (only used if specify_adapter = True).

        -threads
        Number of CPU cores used by fastp.

    Tool: STAR (-tools.STAR)

        -genomeLoad
        How STAR loads genome into memory.
            Load = standard/safe option (use for HPC)
            LoadAndKeep = keeps genome in RAM between runs (only use on same node)

        -outFilterMultimapNmax
        Maximum allowed number of multi-mapping locations. Reads exceeding this are discarded.

        -outFilterMismatchNoverLmax
        max mismatch score before read is discarded (optional, we use star default usually) just leave it blank unless you have a good reason
        0.2 would mean that the alignment must be at least 80% correct

        -outSAMtype
        Output file format and sorting method. Examples:
            BAM Unsorted
            BAM SortedByCoordinate
            BAM SortedByName
            SAM Unsorted

        -twopassMode
            Basic = improves splice junction detection
            None = faster but less accurate

        -file_type
        Command used for reading FASTQs.
        Typically "zcat" for .fastq.gz files or "cat" for uncompressed FASTQ.

        -sjdbOverhang
        Should be read_length - 1.
            you can test your read lengths with this command in linux terminal (WSL):
                zcat /path/to/file.fastq.gz | awk 'NR%4==2 {print length(0$)}' | head -100 | sort | unique -c'
            which prints out paris of numbers, the first is the number of reads that have length equal to the second number
            this looks at the first 100 reads, you can replace head -100 with shuf -n 100 to do 100 random reads from the file but this take a LOT longer to process

        -alignIntronMax
        Maximum allowed intron size. 1,000,000 recommended for mammalian genomes.

        -outReadsUnmapped
        Controls storage of unmapped reads.
            False = do not save
            Within = store inside BAM
            Fastx = save to separate FASTQ output

        -threads
        Number of CPU cores used during STAR alignment.

    Tool: samtools (-tools.samtools)

        -sortMemory
        Maximum memory per thread for samtools sort.
        Total RAM used = sortMemory × threads.

        -minMapQuality
        Minimum MAPQ score for filtering.
            40 = extremely strict, used in variant calling not RNA-seq
            30 = require ~99.9% confidence in map quality
            20 = requrie ~99% confidence in map quality (default)
            10 = requires 90% confidence in map quality
            0 = accepts anyting
            typically 20 or 30 is used, by default we will use 20

        -filter1 and -filter2
        Bitwise flags for filtering reads:
            0x2 = remove improperly paired
            0x4 = remove unmapped
            0x400 = remove duplicates
            0x404 = remove unmapped and duplicates
            Leave blank or False to skip filtering.

        -chromosome
        Restrict samtools operations to a specific chromosome by number (optional).
        Example: 2 

        -chromosomeRange
        Coordinate range on that chromosome (requires chromosome to be set).
        Example: 100000-200000.

        -threads
        Number of CPU cores used for sorting, filtering, and indexing.

    Tool: featureCounts (-tools.featureCounts)

        -gtf_attr_type
        Attribute in GTF used to assign read counts. Usually gene_id.

        -feature_type
        Genomic feature type to count (exon, gene, CDS).

        -strand_specific
        Library strandedness:
            0 = unstranded
            1 = forward-stranded
            2 = reverse-stranded (default, standard for modern Illumina RNA-seq)

        -isPairedEnd
            True = paired-end sequencing
            False = single-end sequencing

        -largestOverlap (optional)
        If True, assign ambiguous reads to the feature with the largest overlap.
        False by default

        -fracOverlap (optional)
        Minimum fraction of a feature that must overlap a read.

        -ignoreDup
        Whether to ignore PCR duplicates.

        -countFraction
        Whether to count multi-mapping reads fractionally.

        -threads
        Number of CPU cores used by featureCounts.

====================================================== Testing =====================================================

Test is split into two parts

    1. FastP/STAR/SamTools