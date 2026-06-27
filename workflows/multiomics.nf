#!/usr/bin/env nextflow

/*
 * Multiomics — Nextflow DSL2 wrapper
 *
 * Runs the multiomics binary inside a container and publishes all outputs
 * (HTML report + JSON) to the directory specified by ``params.output``.
 *
 * Usage
 * -----
 *   nextflow run multiomics.nf \
 *     --genomics     sample.vcf.gz \
 *     --transcriptomics counts.tsv \
 *     --epigenomics  methylation.bed \
 *     --output       results/
 *
 * Optional parameters
 * -------------------
 *   --atac       peaks.narrowPeak
 *   --cnv        cnv.vcf
 *   --gmt        pathways.gmt
 *   --threads    8
 *   --raw_counts (flag — pass as --raw_counts true)
 *   --json_only  (flag — emit JSON only, no HTML/TUI)
 *   --no_ml      (flag — skip PCA and correlation matrix)
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Default parameters — override from the command line or a params file.
// ---------------------------------------------------------------------------

params.genomics         = null          // required: VCF file
params.transcriptomics  = null          // required: expression matrix TSV
params.epigenomics      = null          // required: BED methylation file
params.atac             = null          // optional: ENCODE narrowPeak BED
params.cnv              = null          // optional: VCF with CN INFO field
params.gmt              = null          // optional: GMT pathway file
params.output           = "multiomics_out"
params.threads          = Runtime.runtime.availableProcessors()
params.raw_counts       = false
params.json_only        = false
params.no_ml            = false

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------

def require_param(name, value) {
    if (value == null) {
        error "ERROR: --${name} is required but was not provided."
    }
}

require_param("genomics",        params.genomics)
require_param("transcriptomics", params.transcriptomics)
require_param("epigenomics",     params.epigenomics)

// ---------------------------------------------------------------------------
// Process: BIOOMICS
// ---------------------------------------------------------------------------

process BIOOMICS {
    tag { genomics_file.simpleName }

    label "multiomics"

    // Publish every output file to params.output.
    publishDir params.output, mode: "copy", overwrite: true

    // Container image — overridden per-profile in nextflow.config.
    container "ghcr.io/diladeniz/multiomics:latest"

    input:
    path genomics_file
    path transcriptomics_file
    path epigenomics_file
    path atac_file
    path cnv_file
    path gmt_file

    output:
    path "multiomics_out/report.html",               emit: html,   optional: true
    path "multiomics_out/multiqc_multiomics.json",     emit: json
    path "multiomics_out/**",                        emit: all_outputs

    script:
    // Build the command incrementally so that optional flags are only
    // included when the corresponding input is supplied.
    def atac_arg       = (atac_file.name != "NO_ATAC")        ? "--atac ${atac_file}"              : ""
    def cnv_arg        = (cnv_file.name  != "NO_CNV")         ? "--cnv ${cnv_file}"                : ""
    def gmt_arg        = (gmt_file.name  != "NO_GMT")         ? "--gmt ${gmt_file}"                : ""
    def raw_arg        = params.raw_counts                     ? "--raw-counts"                     : ""
    def json_arg       = params.json_only                      ? "--json"                           : ""
    def no_ml_arg      = params.no_ml                         ? "--no-ml"                          : ""

    """
    multiomics \\
        --genomics        ${genomics_file} \\
        --transcriptomics ${transcriptomics_file} \\
        --epigenomics     ${epigenomics_file} \\
        ${atac_arg} \\
        ${cnv_arg} \\
        ${gmt_arg} \\
        --output          multiomics_out \\
        --threads         ${params.threads} \\
        ${raw_arg} \\
        ${json_arg} \\
        ${no_ml_arg}
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    // Required inputs — each becomes a single-element channel.
    ch_genomics        = Channel.fromPath(params.genomics,        checkIfExists: true)
    ch_transcriptomics = Channel.fromPath(params.transcriptomics, checkIfExists: true)
    ch_epigenomics     = Channel.fromPath(params.epigenomics,     checkIfExists: true)

    // Optional inputs — use a sentinel "NO_*" file when not provided so
    // that Nextflow can still resolve the input declaration.
    ch_atac = params.atac != null
        ? Channel.fromPath(params.atac, checkIfExists: true)
        : Channel.value(file("NO_ATAC"))

    ch_cnv = params.cnv != null
        ? Channel.fromPath(params.cnv, checkIfExists: true)
        : Channel.value(file("NO_CNV"))

    ch_gmt = params.gmt != null
        ? Channel.fromPath(params.gmt, checkIfExists: true)
        : Channel.value(file("NO_GMT"))

    BIOOMICS(
        ch_genomics,
        ch_transcriptomics,
        ch_epigenomics,
        ch_atac,
        ch_cnv,
        ch_gmt,
    )

    // Emit output paths for downstream workflows to consume.
    BIOOMICS.out.json.view { json -> "Multiomics JSON report: ${json}" }
    BIOOMICS.out.html.view { html -> "Multiomics HTML report: ${html}" }
}
