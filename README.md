# Nanopore directRNAsequencing analysis pipeline
Mapping Nanopore directRNAsequencing data

Prerequisites:
1. Java 11 or higher
2. Miniconda
3. Nextflow 
4. Docker

### Note
- Modify the 'nextflow.config' file for Linux systems by changing the executor to 'local' and setting the appropriate paths for your environment. For example:
- Modify the 'nextflow.config' file for linux system

### How to run?
1. Install Miniconda specific to your OS and then install NextFlow
2. Install Docker specific to your OS
3. Download the requisite Docker images using the following command
```
$ docker pull community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797
$ docker pull community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c
```
4. Now cd to `Nanopore_workflow`
5. Then run the command 'nextflow run nanopore.nf'


# Resources
- Read Nextflow instructions [here](https://www.nextflow.io/).
- Test fastq data taken from GitHub repo [here](https://github.com/novoalab/Best_Practices_dRNAseq_analysis/blob/master/README.md).
- Test fastq data taken from [GitHub repo](https://github.com/novoalab/Best_Practices_dRNAseq_analysis/blob/master/README.md)
- Docker containers used from https://seqera.io/containers/
- Biocontainers available at [Quay.io](https://quay.io/).


## Citation
Begik O*, Lucas MC*, Ramirez JM, Milenkovic I, Cruciani S, Vieira HGS, Medina R, Liu H, Sas-Chen A, Mattick JS, Schwartz S and Novoa EM. Decoding ribosomal RNA modification dynamics at single molecule resolution. bioRxiv 2020. doi: https://doi.org/10.1101/2020.07.06.189969



