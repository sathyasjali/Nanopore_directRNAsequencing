# Nanopore directRNAsequencing analysis pipeline
## Overview
This Nextflow pipeline processes sequencing data through various quality control, alignment, variant calling, and visualization steps. It utilizes several bioinformatics tools and modules to ensure comprehensive data processing.

## Workflow Steps
1. **Quality Control (FASTQC, MULTIQC)**: Assesses raw sequencing reads.
2. **Read Filtering (NANOFILT)**: Filters and trims sequencing reads.
3. **Secondary Quality Control (FASTQC)**: Runs after read filtering.
4. **Alignment (MINIMAP2, SAMTOOLS, EMBOSS Needle)**: Maps filtered reads to the reference genome and processes alignment files.
5. **Variant Calling (BCFTOOLS)**: Identifies genetic variants.
6. **Coverage Analysis (MOSDEPTH)**: Computes genome coverage statistics.
7. **Data Visualization (PLOTTING, READLENGTH_HISTOGRAM, ERROR_PLOT)**: Generates plots.
8. **Mapping Statistics Analysis (MAPPING_STATS_ANALYSIS, MERGE_CSV_FILES)**: Computes and consolidates mapping statistics.


## Prerequisites

1. **Java 11 or higher**
2. **Miniconda**
3. **Nextflow**
4. **Docker**


## Requirements
### Dependencies
Ensure that the following software and dependencies are installed:
- [Nextflow](https://www.nextflow.io/)
- Docker or Singularity for containerized execution
- Required bioinformatics tools available via containers:

| **Tool**       | **Version**  | **Container Source** | **Learn More** |
|---------------|------------|--------------------|----------------|
| **FASTQC**     | 0.11.9  | `quay.io/biocontainers/fastqc:0.11.9--0` | [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| **MULTIQC**    | 1.9     | `quay.io/biocontainers/multiqc:1.9--py_1` | [MULTIQC](https://multiqc.info/) |
| **NANOFILT**   | Latest  | `community.wave.seqera.io/library/bcftools_fastqc_minimap2_samtools_pruned:07a664876843239a` | [NANOFILT](https://github.com/wdecoster/nanofilt) |
| **MINIMAP2**   | Latest  | `community.wave.seqera.io/library/bcftools_emboss_minimap2_samtools:09b81fb0810b8e21` | [MINIMAP2](https://github.com/lh3/minimap2) |
| **SAMTOOLS**   | Latest  | `community.wave.seqera.io/library/bcftools_emboss_minimap2_samtools:09b81fb0810b8e21` | [SAMTOOLS](http://www.htslib.org/) |
| **BCFTOOLS**   | Latest  | `community.wave.seqera.io/library/bcftools_emboss_minimap2_samtools:09b81fb0810b8e21` | [BCFTOOLS](http://www.htslib.org/) |
| **EMBOSS**     | Latest  | `community.wave.seqera.io/library/bcftools_emboss_minimap2_samtools:09b81fb0810b8e21` | [EMBOSS](http://emboss.sourceforge.net/) |
| **MOSDEPTH**   | 0.3.10  | `community.wave.seqera.io/library/mosdepth:0.3.10--259732f342cfce27` | [MOSDEPTH](https://github.com/brentp/mosdepth) |

### Input Data
- **FASTQ Files**: The pipeline expects input FASTQ sequencing files located at:
  ```
  data/reads/*.fastq
  ```
- **Reference Genome**: The reference genome must be placed in:
  ```
  data/reference/*.fa
  ```

### Configuration
Modify the `params` section in the Nextflow script to adjust input file locations and execution parameters.


### Note
- Modify the `nextflow.config` file for Linux systems by changing the executor to 'local' and setting the appropriate paths for your environment. For example:
- Modify the `nextflow.config` file for linux system

### How to run?
1. Install Miniconda specific to your OS and then install NextFlow
2. Install Docker specific to your OS
3. Now `cd` to `Nanopore_directRNAsequencing`
4. Then run the command `nextflow run nanopore.nf`
5.  Nextflow skips re-executing completed steps and resumes execution from the point of failure or where changes occurred.`nextflow run nanopore.nf -resume`

## Output Files
The pipeline generates multiple output files in the `results` directory:
- **Quality Control Reports (FASTQC)** → `results/fastqc/`
- **Filtered Reads (NANOFILT)** → `results/nanofilt/`
- **Alignment Files (MINIMAP2, SAMTOOLS)** → `results/minimap2/`
- **Variant Calls (BCFTOOLS)** → `results/bcftools/`
- **Coverage Statistics (MOSDEPTH)** → `results/mosdepth/`
- **Plots & Graphs (PLOTTING, HISTOGRAM, ERROR_PLOT)** → `results/plots/`
- **Final Merged CSV Report (MAPPING_STATS_ANALYSIS)** → `results/mapping_stats.csv`
- **HTML Report** → `results/report_summary.html`

## Generate Report
The pipeline includes a final step to generate a report summarizing all key outputs. To generate the report, run:
```sh
nextflow run nanopore.nf -profile docker --generate_report true
```
This will output a structured report file at `results/report_summary.html`, containing:
- Summary of quality control metrics
- Read filtering statistics
- Alignment and variant calling results
- Coverage analysis overview
- Key plots and visualizations

## Docker Configuration
The pipeline uses Docker containers for reproducibility. If running with Docker, ensure Docker is installed and running. The pipeline will automatically pull necessary images.

To manually build the required environment, use:
```sh
docker pull quay.io/biocontainers/fastqc:0.11.9--0
```

### Folder structure ##

```
Nextflow Workflow
├── Quality Control
│   ├── FASTQC (Raw Reads)
│   ├── MULTIQC (FASTQC Reports)
│   └── FASTQC (Filtered Reads)
├── Read Filtering
│   ├── NANOFILT
│   └── FASTQC (Post-Filtering)
├── Alignment
│   ├── MINIMAP2 (Filtered Reads → Reference Genome)
│   ├── SAMTOOLS (Sort, Index, Process BAM)
│   └── EMBOSS Needle (Alignment Refinement)
├── Variant Calling
│   ├── BCFTOOLS (Variant Identification)
│   ├── Extract VCF Output
│   └── ERROR_PLOT (Error Visualization)
├── Coverage Analysis
│   ├── MOSDEPTH (Depth Analysis)
│   └── Coverage Extraction
├── Data Visualization
│   ├── PLOTTING (General Plots)
│   ├── READLENGTH_HISTOGRAM (Read Length Analysis)
│   ├── ERROR_PLOT (Error Visualization)
│   └── MULTIQC (Final Report)
├── Statistics Analysis
│   ├── MAPPING_STATS_ANALYSIS (Mapping Statistics)
│   ├── MERGE_CSV_FILES (Consolidate Stats)
│   └── Final Report Generation
└── Output Storage
    ├── FASTQC Reports → results/fastqc/
    ├── Filtered Reads → results/nanofilt/
    ├── Alignment Files → results/minimap2/
    ├── Variant Calls (VCF) → results/bcftools/
    ├── Coverage Statistics → results/mosdepth/
    ├── Plots & Graphs → results/plots/
    └── Final Merged CSV Report → results/mapping_stats.csv

```

## Troubleshooting
### Common Issues & Fixes
1. **No FASTQ files found**
   - Ensure input files are placed in `data/reads/`
   - Verify file paths in `params.fastq_files`

2. **Reference genome not found**
   - Ensure reference genome files are in `data/reference/`
   - Check `params.reference_dir`

3. **Missing dependencies**
   - Use `nextflow info` to verify module availability
   - Ensure required containers are accessible


# Resources
- Read Nextflow instructions [here](https://www.nextflow.io/).
- Test fastq data taken from GitHub repo [here](https://github.com/novoalab/Best_Practices_dRNAseq_analysis/blob/master/README.md).
- Test fastq data taken from [GitHub repo](https://github.com/novoalab/Best_Practices_dRNAseq_analysis/blob/master/README.md)
- Docker containers used from https://seqera.io/containers/
- Biocontainers available at [Quay.io](https://quay.io/).


## Citation
Begik O*, Lucas MC*, Ramirez JM, Milenkovic I, Cruciani S, Vieira HGS, Medina R, Liu H, Sas-Chen A, Mattick JS, Schwartz S and Novoa EM. Decoding ribosomal RNA modification dynamics at single molecule resolution. bioRxiv 2020. doi: https://doi.org/10.1101/2020.07.06.189969
