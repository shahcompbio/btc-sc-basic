#!/usr/bin/env nextflow

log.info """\
BTC - SINGLE-CELL BASIC PIPELINE
===================================
Project parameters:
- Project Name          : ${params.project_name}
- Sample CSV            : ${params.sample_csv}
- Meta-data             : ${params.meta_data}
- Genome                : ${params.genome}
- Fasta                 : ${params.fasta}
- Annotation            : ${params.annotation}
"""

/*
    Download Human Indexes
*/

process DOWNLOAD_HG_INDEX {
    /* 
    Downloading genome indexes from BTC Bucket
    https://storage.googleapis.com/btc-dshub-pipelines/scRNA/refData/GRCh38.tar 
    */
    
    label "Download_Human_Indexes"

    input:
        val(genome) 

    output:
        path("indexes/${genome}"), emit: index

    script:
        indexes = ["GRCh38" : "https://storage.googleapis.com/btc-dshub-pipelines/scRNA/refData/GRCh38.tar"]

        """
        wget -O ${genome}.tar ${indexes[genome]}
        mkdir ./indexes
        tar -xvf ${genome}.tar -C ./indexes
        rm -rf ${genome}.tar
        """
}

/*
    Custom reference (Run more tests)
*/

process CUSTOM_REFERENCE {
    /* Description */

    container 'oandrefonseca/scaligners:1.0'
    label "Custom_Reference" 

    input:
        val(genome)
        path(fasta)
        path(annotation)

    output:
        path("${annotation.baseName}.filtered.gtf")
        path("indexes/${genome}")

    script:
        anno_filtered = "${annotation.baseName}.filtered.gtf"
        """
        cellranger \\
            mkgtf \\
            ${annotation} \\
            ${anno_filtered} \\
            --attribute=gene_type:protein_coding

        cellranger \\
            mkref \\
                --genome="indexes/${genome}" \\
                --fasta=${fasta} \\
                --genes=${anno_filtered} \\
                --nthreads=${params.cpus}
        """
}

/*
    Sample Alignment
*/

process SAMPLE_ALIGNMENT {
    /* Description */

    container 'oandrefonseca/scaligners:1.0'
    label "Sample_Alignment"

    tag "Processing ${sample_id}"
    publishDir "${params.project_name}/", mode: 'copyNoFollow'

    input:
        tuple val(sample_id), val(prefix), path(first_read), path(second_read)
        path(reference)

    output:
        tuple val(sample_id), path("${sample_id}/outs/*"), emit: cell_out

    script:
        """
        cellranger \\
            count \\
            --id='${sample_id}' \\
            --fastqs=. \\
            --transcriptome=${reference} \\
            --include-introns=false \\
            --no-bam \\
            --sample=${prefix} \\
            --localcores=${params.cpus} \\
            --localmem=${params.memory}
        """

}

/*
    Sample and Cell quality control
*/

process SAMPLE_CELL_QC {
    /* Description */

    container 'oandrefonseca/scpackages:1.1'
    label "Sample_and_Cell_QC"

    tag "QC ${sample_id}"
    publishDir "${params.project_name}/${sample_id}", mode: 'copyNoFollow'

    input:
        tuple val(sample_id), path(matrices), path(csv_metrics), path(meta_data)
        path(scqc_script)

    output:
        tuple val(sample_id), path("objects/*"), path("log/*.txt"), emit: status
        path("${sample_id}_metrics_upgrade*.csv"), emit: metrics
        path("${sample_id}_report.html")
        path("figures/*")

    script:
        """
        #!/usr/bin/env Rscript

        # Getting run work directory
        here <- getwd()

        # Rendering Rmarkdown script
        rmarkdown::render("${scqc_script}",
            params = list(
                project_name = "${params.project_name}",
                sample_name = "${sample_id}",
                meta_data = "${meta_data}",
                matrices = "${matrices}",
                csv_metrics = "${csv_metrics}",
                thr_estimate_n_cells = ${params.thr_estimate_n_cells},
                thr_mean_reads_per_cells = ${params.thr_mean_reads_per_cells},
                thr_median_genes_per_cell = ${params.thr_median_genes_per_cell},
                thr_median_umi_per_cell = ${params.thr_median_umi_per_cell},
                thr_nFeature_RNA_min = ${params.thr_nFeature_RNA_min},
                thr_nFeature_RNA_max = ${params.thr_nFeature_RNA_max},
                thr_percent_mito = ${params.thr_percent_mito},
                thr_n_observed_cells = ${params.thr_n_observed_cells},
                workdir = here
            ), 
            output_dir = here,
            output_file = "${sample_id}_report.html")           
        
        """
}

/*
    Merging count matrices
*/

process QUALITY_TABLE {
    /* Description */

    container 'oandrefonseca/scpackages:1.1'
    label "Quality_table"
    
    publishDir "${params.project_name}", mode: 'copyNoFollow'

    input:
        path(project_metrics)
        path(qc_table_script)

    output:
        path("${params.project_name}_project_metric_report.html")
    
    script:
        """
        #!/usr/bin/env Rscript

        # Getting run work directory
        here <- getwd()

        # Rendering Rmarkdown script
        rmarkdown::render("${qc_table_script}",
            params = list(
                project_name = "${params.project_name}",
                input_metrics_report = "${project_metrics.join(';')}",
                workdir = here
            ), 
            output_dir = here,
            output_file = "${params.project_name}_project_metric_report.html")
        """        

}

/*
    Merging count matrices
*/

process MERGE_AND_NORMALIZE {
    /* Description */

    container 'oandrefonseca/scpackages:1.1' 
    label "Merge_and_Normalize"

    publishDir "${params.project_name}", mode: 'copyNoFollow'

    input:
        path(good_qc_matrices)
        path(merge_script)

    output:
        path("${params.project_name}_normalize_object.RDS"), emit: project_rds
        path("${params.project_name}_normalize_report.html")
        path("figures/*")

    script:
        """
        #!/usr/bin/env Rscript

        # Getting run work directory
        here <- getwd()

        # Rendering Rmarkdown script
        rmarkdown::render("${merge_script}",
            params = list(
                project_name = "${params.project_name}",
                input_qc_approved = "${good_qc_matrices.join(';')}",
                workdir = here
            ), 
            output_dir = here,
            output_file = "${params.project_name}_normalize_report.html")
        """
}

/*
    Batch correction
*/

process BATCH_CORRECTION {

    container 'oandrefonseca/scpackages:1.1' 
    label "Batch_Correction"

    publishDir "${params.project_name}", mode: 'copyNoFollow'

    input:
        path(project_object)
        path(batch_script)

    output:
        path("${params.project_name}_batch_object_*.RDS"), emit: project_rds
        path("${params.project_name}_batch_report.html")
        path("figures/*")

    script:
        """
        #!/usr/bin/env Rscript

        # Getting run work directory
        here <- getwd()

        # Rendering Rmarkdown script
        rmarkdown::render("${batch_script}",
            params = list(
                project_name = "${params.project_name}",
                project_object = "${project_object}",
                input_target_variables = 'batch',
                workdir = here,
                timestamp = "${workflow.runName}"
            ), 
            output_dir = here,
            output_file = "${params.project_name}_batch_report.html")           

        """
}

/*
    Cell clustering

*/

process CELL_CLUSTERING {
  
    container 'oandrefonseca/scpackages:1.1' 
    label "Cell_clustering"

    publishDir "${params.project_name}", mode: 'copyNoFollow'

    input:
        path(project_object)
        path(cluster_script)

    output:
        path("${params.project_name}_cluster_object_*.RDS"), emit: project_rds
        path("${params.project_name}_cluster_report.html")
        path("figures/*")
        path("data/*")

    script:
        """
        #!/usr/bin/env Rscript

        # Getting run work directory
        here <- getwd()

        # Rendering Rmarkdown script
        rmarkdown::render("${cluster_script}",
            params = list(
                project_name = "${params.project_name}",
                project_object = "${project_object}",
                workdir = here,
                timestamp = "${workflow.runName}"

            ), 
            output_dir = here,
            output_file = "${params.project_name}_cluster_report.html")           

        """

}

/*
    Single-cell RNA-Seq pipeline
    ----------------------------
*/

workflow {

    /* Main workflow body */

    if (!params.project_name) { exit 1, 'Input project name not specified!' }
    if (!params.sample_csv) { exit 1, 'Input samples file not specified!' }
    if (!params.meta_data) { exit 1, 'Input meta-data not specified!' }

    if (params.genome == 'GRCh38') {
        DOWNLOAD_HG_INDEX(
            params.genome
        )
    }

    if (params.genome != 'GRCh38') {

        if(params.fasta && params.annotation) {

            CUSTOM_REFERENCE(
                params.genome,
                params.fasta,
                params.annotation
            )

        }

    }
    
    // Rmarkdown scripts 
    scqc_script = "${workflow.projectDir}/bin/01_quality_control.Rmd"
    qc_table_script = "${workflow.projectDir}/bin/02_quality_table_report.Rmd"
    merge_script = "${workflow.projectDir}/bin/03_merge_and_normalize.Rmd"
    batch_script = "${workflow.projectDir}/bin/04_batch_correction.Rmd"
    cluster_script = "${workflow.projectDir}/bin/05_cell_clustering.Rmd"

    // Samples channel
    samples_channel = Channel.fromPath(params.sample_csv)
                             .splitCsv(header:true)

    // Meta-data path
    meta_channel = Channel.fromPath(params.meta_data)

    // Cellranger alignment
    SAMPLE_ALIGNMENT(samples_channel, DOWNLOAD_HG_INDEX.out.index)

    cell_matrices_channel = SAMPLE_ALIGNMENT.out.cell_out
    .map{sample, outs -> [sample, outs.findAll {
        it.toString().endsWith("metrics_summary.csv") || it.toString().endsWith("filtered_feature_bc_matrix")
    }]}
    .map{sample, files -> [sample, files[0], files[1]]}

    cell_matrices_channel = cell_matrices_channel
    .combine(meta_channel)

    // Performing QC steps
    SAMPLE_CELL_QC(cell_matrices_channel, scqc_script)

    // Writing QC check
    quality_report_channel = SAMPLE_CELL_QC.out.metrics
    .collect()
    
    QUALITY_TABLE(quality_report_channel, qc_table_script)

    // Filter poor quality samples
    qc_approved_channel = SAMPLE_CELL_QC.out.status
    .filter{sample, object, status -> status.toString().endsWith('SUCCESS.txt')}
    .map{sample, object, status -> object}
    .collect()

    qc_approved_channel
    .ifEmpty{error 'No samples matched QC expectations.'}
    .view{'Done'}

    // Merging datasets
    MERGE_AND_NORMALIZE(qc_approved_channel, merge_script) 
    
    // Batch Correction and Cell Clustering
    if(params.skip_batch) {

        CELL_CLUSTERING(MERGE_AND_NORMALIZE.out.project_rds, cluster_script)

    } else {

        BATCH_CORRECTION(MERGE_AND_NORMALIZE.out.project_rds, batch_script)
        CELL_CLUSTERING(BATCH_CORRECTION.out.project_rds, cluster_script)

    }

}

workflow.onComplete {

    log.info(workflow.success ? "May the Force be with you!" : "Please check your inputs.")

}

