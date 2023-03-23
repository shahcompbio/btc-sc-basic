# btc-sc-basic

> **What is a minimum viable product?** A minimum viable product (MVP) is a version of a product with just enough features to be usable by early customers who can then provide feedback for future product development. [1](https://en.wikipedia.org/wiki/Minimum_viable_product)

> **What does the MVP version do?** The pipeline was designed to handle sample processing and quality control. It also performs “basic” analysis related to batch correction, clustering, and differential expression. As a result, multiple outputs are generated in each step, including sample-level reports and project-related Seurat objects.


## 1. Installing Nextflow and third-party software

Nextflow can be used on any POSIX-compatible system (Linux, OS X, WSL). It requires Bash 3.2 (or later) and Java 11 (or later, up to 18) to be installed.

```bash

wget -qO- https://get.nextflow.io | bash

```

## 1.1. Docker and Singularity

This tutorial was tested using Docker v20.10.22 and Singularity v3.7.0. Please make sure to have both software available. The `btc-sc-basic` pipeline will require docker images for running.

## 2. Cloning `btc-sc-basic`

```bash

git clone git@github.com:WangLab-ComputationalBiology/btc-sc-basic.git

``` 

## 3. Building containers from scratch

A brief disclaimer regarding non-academic organizations: the DockerHub license is quite expensive. To avoid extra expenses and legal issues, we avoid pulling images from any repository. Currently, the BTC infrastructure is working in alternatives for DockerHub. Meanwhile, let's build the images from scratch.

### 3.1 Running with Docker images

Two Dockerfiles are available in the folder `container`. Each file corresponds to a single Docker image. By running the following commands, we can build Docker images from these recipes. Note the `scAligners.Dockerfile` requires a URL derived from 10x Genomics website. For more explanation, please get in touch.

```bash

docker build -t scaligners:1.0 -f scAligners.Dockerfile .

``` 

```bash

docker build -t scpackages:1.1 -f scPackages.Dockerfile .

``` 

### 3.2 Running with Singularity images

Docker solutions are "supposedly" dangerous for HPC environments. Therefore, we must use Singularity (v3.7.0) to run the single-cell pipeline on computational clusters. For this purpose, we suggest downloading the pre-built images from our BTC bucket [**UNDER REQUEST**].

Additionally, to optimize the Nextflow performance on clusters, we can export `NXF_SINGULARITY_CACHEDIR` to shared directory. This procedure will share Singularity images across pipelines/runs. It could be done the following way:

```bash

export NXF_SINGULARITY_CACHEDIR=$SCRATCH/singularity

``` 

**This instruction should be written on** `.bash_profile`.

---

## 4.  Running single-cell pipeline

### 4.1. Preparing inputs

The pipeline requires two inputs, a sample table and meta-data. Both files should follow a **mandatory** format as described below.

1. Sample table should be a `csv` with four columns, sample, prefix, fastq_{1,2}. The column `sample` will be associated with all reports across the pipeline. Also, it will be required to combine the meta-data on the Seurat object. Prefix and fastq columns are mandatory for `Cellranger`.

|   sample_id  | prefix |  fastq_1 | fastq_2 |
|:---------:|:------:|:--------:|:-------:|
| BTC-AA | BTC-AA-001 | path/to/BTC-AA-001_S1_L001_R1_001.fastq.gz | path/to/BTC-AA-001_S1_L001_R2_001.fastq.gz |
| BTC-AB | BTC-AB-001 | path/to/BTC-AB-001_S1_L001_R1_001.fastq.gz | path/to/BTC-AB-001_S1_L001_R2_001.fastq.gz |
| BTC-AC | BTC-AC-001 | path/to/BTC-AC-001_S1_L001_R1_001.fastq.gz | path/to/BTC-AC-001_S1_L001_R2_001.fastq.gz |

**Important, note the `prefix` column has to match fastq files**

2. Meta-data (`.csv`) should contain variables related to experimental design (e.g., batch and cell sorting status). Additionally  it could have more biological information about that sample. The batch variable will be used to correct the effect. In this pipeline version, we are executing a correction based on a single variable. Finally, the sort column will be used in future editions as an additional step for the QC process.

| sample_id | batch | sort | ... |
|:------:|:-----:|:----:|:---:|
| BTC-AA | CORE_A | CD45+ | ... |
| BTC-AB | CORE_A | CD45+ | ... |
| BTC-AC | CORE_B | CD45- | ... |

### 4.2. Preparing parameters for pipeline execution

To run the pipeline, we must input three parameters: `project_name`, `sample_csv`, and `meta_data`. In addition, for those using the HPC environment, we need to set up a **profile**. Briefly, profiles are specific settings related to the `job scheduler` on the HPC. For instance, each HPC can rely on different engines to schedule jobs, such as SLURM, TORQUE, and LSF. Additionally, profiles can be used to store parameters, i.e., for testing purposes.

```bash

nextflow run single_cell_basic.nf --project_name <PROJECT> --sample_csv <path/to/sample_table.csv> --meta_data <path/to/meta_data.csv> -resume -profile <PROFILE>

```

**Note**:  There is a semantic difference between the `nextflow` (-) and the `pipeline` (--) parameters. For instance, double-dashed parameters are restricted to pipeline logic, e.g., coded to change filters and thresholds on the single-cell analysis. Finally, the pipeline will generate a folder based on `--project_name`, which will hold our results after the end of the execution. The `-resume` command will leverage the Nextflow caching system to provide re-entrance.

### 4.3. An example of a SLURM profile

The profiles should be written on the `nextflow.config` file. For HPCs, loading the **singularity/3.7.0** will be mandatory. We can add information about queue size, memory, and cpus. 

```java

slurm {

    module = 'singularity/3.7.0'

    singularity {
        enabled = true
    }

    process {
        executor = 'slurm'
        cpu = 24
        queue = 'medium'
        memory = '64 GB'
    }

    params {
        cpus = 24
        memory = 64
    }

}

```

### 4.3. parameters per process/step

Currently, the MVP version has several parameters per process/analysis/data flow. Please check out, the parameters list and default values below.

### Alignment-related parameters
* `genome` = 'GRCh38' # Pre-computed index based on basic gencode v43 (Default)
* `fasta` = false # OPTIONAL: Path to a fasta file
* `annotation` = false # OPTIONAL: Path to a GTF file

### Quality control-related parameters
* `thr_estimate_n_cells` = 300
* `thr_mean_reads_per_cells` = 25000
* `thr_median_genes_per_cell` = 900
* `thr_median_umi_per_cell` = 1000
* `thr_nFeature_RNA_min` = 300
* `thr_nFeature_RNA_max` = 7500
* `thr_percent_mito` = 25
* `thr_n_observed_cells` = 300

### Dimensionality reduction and batch correction
* `input_target_variables` = 'batch' # The Harmony model will be performed using it as the target variable. It can be multiple columns on your meta-data.
* `thr_npc` = 'auto' # Adjusted based on the number of cells. The user can input a integer threshold.

### Clustering- and UMAP-related parameters
* `input_features_plot` = 'LYZ;CCL5;IL32;PTPRCAP;FCGR3A;PF4;PTPRC' # A list of markers to be plotted
* `input_group_plot`: 'batch' # Generates UMAP grouped by one or more variables. 
* `run_deg` = FALSE
* `thr_quantile` = 'q01'
* `thr_resolution` = 0.25 # A single resolution value per execution

### 4.4. Skipping analytical steps (processes) on the pipeline

In some situations, the user might want to skip processes. Currently, users can avoid batch correction analysis by adding the `--skip_batch`. That might help analyze the impact of the batch correction on the final clustering (before and after batch correction).

### 4.5. Shorten command-line 

A huge command line can be painful. Luckily, we can shorten instructions by using Nextflow `-params-file`. It will be a JSON file containing all parameters for that specific run. In the case of testing parameters, it could be a good practice to generate multiple files (e.g., PARAMS_TEST_01.json, PARAMS_TEST_02.json).

```json

{
 "project_name": "BTC-DISEASE-00",
 "sample_csv": "path/to/sample_table.csv",
 "meta_data": "path/to/meta_data.csv"
}

```

In addition, we could add custom QC parameters:

```json

{
 "project_name": "BTC-DISEASE-00",
 "sample_csv": "path/to/sample_table.csv",
 "meta_data": "path/to/meta_data.csv",
 "thr_mean_reads_per_cells": 10000
}

```

```bash

nextflow run single_cell_basic.nf -params-file <PARAMS.json> -resume -profile <PROFILE>

```

---

## 5. Expected outputs

In a successful execution, the `btc-sc-basic` should create a folder structure similar to the one described below. The files with TIMESTAMP are created based on the nextflow RUN NAME [1](https://www.nextflow.io/docs/latest/tracing.html). Theoretically, it should improve traceability across executions by levering Nextflow properties.

```

PROJECT_NAME
├── SAMPLE_1
│   ├── figures
│   ├── log
│   │   └── SUCCESS.txt
│   ├── objects
│   └── outs
│       ├── ...
│       ├── filtered_feature_bc_matrix
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   └── matrix.mtx.gz
│       └── raw_feature_bc_matrix
├── SAMPLE_N
│   ├── figures
│   ├── log
│   │   └── FIXABLE.txt
│   ├── objects
│   └── outs
│       ├── ...
│       ├── filtered_feature_bc_matrix
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   └── matrix.mtx.gz
│       └── raw_feature_bc_matrix
├── data
│   └── deg
│   │   ├── PROJECT_NAME_deg_analysis_<TIMESTAMP>.RDS
├── figures
│   ├── batch
│   ├── clustered
│   │   ├── UMAP_FEATURED_<TIMESTAMP>.pdf
│   │   ├── UMAP_GROUPED_1_<TIMESTAMP>.pdf
│   │   ├── UMAP_GROUPED_2_<TIMESTAMP>.pdf
│   │   └── UMAP_MAIN_<TIMESTAMP>.pdf
│   └── normalized
│       └── Elbow_plots_<TIMESTAMP>.pdf
├── PROJECT_NAME_batch_object_<TIMESTAMP>.RDS
├── PROJECT_NAME_batch_report.html
├── PROJECT_NAME_clustering_report.html
├── PROJECT_NAME_clustering_object_<TIMESTAMP>.RDS
├── PROJECT_NAME_normalize_report.html
└── PROJECT_NAME_normalize_object.RDS

```