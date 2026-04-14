// ================================================================
// modules/multivariate_analysis.nf — CDC 1.0.0
// Multivariate statistical analysis module:
//   Phase 1 · Hypothesis testing (MANOVA proxy, Hotelling T²)
//   Phase 2 · Regression (OLS, Ridge, Lasso, EN, PCR, PLS, VIF)
//   Phase 3 · Multiple-testing correction (Bonferroni, Holm, BH, BY, Storey)
//   Phase 4 · Causal inference (2SLS/IVW, MR-Egger, Cochran Q, LOO)
// Author : Nadeem Khan, INRS-CAFSB
// ================================================================

process MULTIVARIATE_ANALYSIS {
    label 'medium'
    tag "multivariate_all_cohorts"

    publishDir "${params.out_dir}/multivariate", mode: 'copy', overwrite: true

    // Environment (conda or container) is set exclusively in nextflow.config

    input:
    path sumstat_files   // all *_sumstats.tsv.gz files staged into work dir
    val  author
    val  institute

    output:
    path "phase1_hypothesis/static/**",    emit: phase1_static,  optional: true
    path "phase2_regression/static/**",    emit: phase2_static,  optional: true
    path "phase3_mtc/static/**",           emit: phase3_static,  optional: true
    path "phase4_causal/static/**",        emit: phase4_static,  optional: true
    path "phase1_hypothesis/*.tsv",        emit: phase1_tables
    path "phase2_regression/*.tsv",        emit: phase2_tables
    path "phase3_mtc/*.tsv",               emit: phase3_tables
    path "phase4_causal/*.tsv",            emit: phase4_tables
    path "phase1_hypothesis/*.xlsx",       emit: phase1_excel,   optional: true
    path "phase2_regression/*.xlsx",       emit: phase2_excel,   optional: true
    path "phase3_mtc/*.xlsx",              emit: phase3_excel,   optional: true
    path "phase4_causal/*.xlsx",           emit: phase4_excel,   optional: true
    path "Multivariate_Report.html",       emit: html_report

    script:
    """
    set -euo pipefail

    echo "───────────────────────────────────────────────"
    echo "  CDC 1.0.0 · Multivariate Analysis"
    echo "  work dir : \${PWD}"
    echo "  n sumstat files staged : \$(ls *_sumstats.tsv.gz 2>/dev/null | wc -l)"
    echo "  Python : \$(python3 --version)"
    echo "───────────────────────────────────────────────"

    # Make modules (static_figures.py, excel_tables.py)
    export PYTHONPATH="${projectDir}/bin:\${PYTHONPATH:-}"

    python3 ${projectDir}/bin/multivariate_analysis.py \\
        --sumstats_pattern "*_sumstats.tsv.gz" \\
        --out_dir           "\${PWD}" \\
        --author            "${author}" \\
        --institute         "${institute}"
    """
}
