#!/usr/bin/env nextflow
// ================================================================
// main.nf — CDC 1.0.0 · Multivariate Statistical Analysis Pipeline
// Extends gwas-pipeline-v3 (post-REGENIE) with:
//   · Multivariate hypothesis testing
//   · Multi-model regression + VIF
//   · Multiple-testing correction (Bonferroni, Holm, BH, BY, Storey)
//   · Causal inference (2SLS/IVW, MR-Egger, Cochran Q, LOO)
//
// Author : Nadeem Khan, INRS-Centre Armand-Frappier Santé-Biotechnologie
// GitHub : github.com/nkhan119
// ================================================================

nextflow.enable.dsl = 2

include { MULTIVARIATE_ANALYSIS } from './modules/multivariate_analysis.nf'

// ── Banner ──────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════╗
║  Multivariate GWAS Pipeline · CDC 1.0.0                        ║
║  Author     : ${params.author.padRight(48)}║
║  Institute  : ${(params.institute ?: "").take(48).padRight(48)}║
╠══════════════════════════════════════════════════════════════════╣
║  sumstats   : ${params.sumstats_dir.padRight(48)}║
║  out_dir    : ${params.out_dir.padRight(48)}║
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()

// ── Workflow ─────────────────────────────────────────────────────
workflow {

    // Collect all harmonised sumstat files from the GWAS output tree
    // Supports both the gwas-pipeline-v3 directory layout and flat dirs
    ch_sumstats = Channel
        .fromPath([
            "${params.sumstats_dir}/**/sumstats/*_sumstats.tsv.gz",
            "${params.sumstats_dir}/*_sumstats.tsv.gz"
        ], glob: true)
        .unique()
        .collect()
        .ifEmpty {
            error """
            ╔══════════════════════════════════════════════════════════════╗
            ║  ERROR: No *_sumstats.tsv.gz files found under:            ║
            ║  ${params.sumstats_dir.padRight(60)}║
            ╚══════════════════════════════════════════════════════════════╝
            Hint: Run gwas-pipeline-v3 first, then set --sumstats_dir to
                  the directory containing gwas/**/sumstats/*.tsv.gz
            """
        }

    MULTIVARIATE_ANALYSIS(
        ch_sumstats,
        params.author,
        params.institute ?: "INRS-Centre Armand-Frappier Santé-Biotechnologie"
    )
}

// ── Completion hook ──────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? "SUCCESS ✓" : "FAILED ✗"
    log.info """
╔══════════════════════════════════════════════════════════════════╗
║  Multivariate Pipeline ${status.padRight(40)}║
║  Duration  : ${workflow.duration.toString().padRight(49)}║
║  HTML      : ${(params.out_dir + '/multivariate/Multivariate_Report.html').take(49).padRight(49)}║
║  Tables    : ${(params.out_dir + '/multivariate/phase*/').padRight(49)}║
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()
}
