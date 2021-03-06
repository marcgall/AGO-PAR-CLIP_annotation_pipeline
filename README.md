Pipeline for annotation of raw reads from AGO-PARCLIP sequence data.

**Requirements**
* Conda (from anaconda or miniconda)
* snakemake (accessable from the conda channel *bioconda*).

Everything else is downloaded as needed by snakemake.

*Mac users:* 
* 'wget' must be available on the command-line.
* PARalyzer requires Java 1.6, which is not currently available through conda for mac. Thus, you must make it available on your own machine AND, to make sure that snakemake doesn't try to download java, go to the envs/PARalyzer.yaml file and remove everything below "-bioconda".
*(If you find a conda package for java 1.6 that is available for both linux and mac, please let me know.)*

**Basic Usage**

The basic form of the pipeline is run by the following command inside the directory: "snakemake --use-conda"
This runs the whole pipeline for all SRR files specified in SRR.csv. To utilize multiple cores on your machine, use the "--cores n" argument, replacing "n" with the number of cores.

The target SRR files can be specified manually in the command-line or entered in the SRR.csv file. Only the first column of the SRR.csv file is read by the pipeline - the rest can be used for notes that aid the user. The current SRR file is supplied as an example.

Output files can be specified by adding the name to the end of the snakemake command.

**Advanced Usage**

Specific settings for cutadapt and bowtie can be changed in the VARIABLES section of the Snakemake file. Note, however, that the downloading of human mature miRs from mirbase is currently hardcoded.

The pipeline is "lazy", meaning that if a current named file in the middle of the pipeline is present, it will skip the preceding steps. Thus any step can be circumvented by supplying the correct file. For example, to supply your own fastq files, make the directory "raw_reads" and place your fastq files there. If not specifying the output file that you want in the snakemake command, be sure to also update the SRR.csv file so that the specific run names correspond between all files.

To keep intermediate files in the pipeline, the 'temp()' tag can be removed in the snakemake rules.

To see what steps will be calculated, run 'snakemake -np'

To generate a graph of the pipeline, run "snakemake --dag | dot -Tpdf > example_graph.pdf"
