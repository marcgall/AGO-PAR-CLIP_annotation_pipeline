Pipeline for annotation of raw reads from AGO-PARCLIP sequence data.

Python and Conda must be installed (through anaconda or miniconda). Everything else is downloaded as needed by snakemake.

The basic form of the pipeline is run by the following command:
$ snakemake --use-conda
This runs the whole pipeline for all SRR files specified in SRR.csv. To utilize multiple cores on your machine, use the "--cores n" argument, replacing "n" with the number of cores.

The target SRR files can be specified manually in the command-line or entered in the SRR.csv file. Only the first column of the SRR.csv file is read by the pipeline - the rest can be used for notes that aid the user. The current SRR file is supplied as an example.

Specific settings for cutadapt and bowtie can be changed in the VARIABLES section of the Snakemake file. Note, however, that the downloading of human mature miRs from mirbase is
currently hardcoded.

At any step, the pipeline is superceded by the presence of correctly named
files, by which downloading or processing new files can be avoided. Similarly,
the 'temp()' tag can be removed in the rules if automatic cleanup is not
desired.
