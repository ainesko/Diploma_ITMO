"""Run the complete mutational signature analysis workflow."""

import signature_extraction.run_musical_denovo as musical
import signature_extraction.run_sigprofiler_denovo as sigprofiler
import signature_extraction.musical_cosmic_assignment as cosmic

from analysis.metrics import (
    compute_similarity_matrix,
    match_signatures,
    compute_explained_variance_table,
)

from analysis.visualization import (
    plot_similarity_heatmap,
    plot_explained_variance,
    plot_cosine_similarity_dotplot,
)


INPUT_MATRIX = "data/example_matrix.tsv"


# MuSiCal de novo signature extraction
musical.run_musical_denovo(
    input_matrix=INPUT_MATRIX,
    output_dir="results/musical/de_novo",
    sig_type="SBS96",
    min_components=1,
    max_components=4,
    n_replicates=20,
    ncpu=1,
)


# COSMIC assignment for MuSiCal signatures
cosmic.run_musical_cosmic_assignment(
    signature_matrix="results/musical/de_novo/musical_W_de_novo.csv",
    output_dir="results/musical/cosmic_assignment",
    sig_type="SBS96",
)


# SigProfilerExtractor de novo signature extraction
sigprofiler.run_sigprofiler(
    input_type="matrix",
    input_data=INPUT_MATRIX,
    output_dir="results/sigprofiler/de_novo",
    reference_genome="GRCh38",
    context_type="96",
    minimum_signatures=1,
    maximum_signatures=4,
    nmf_replicates=100,
    cpu=1,
)


# Compare MuSiCal and SigProfiler de novo signatures
similarity_matrix = compute_similarity_matrix(
    signatures_a="results/musical/de_novo/musical_W_de_novo.csv",
    signatures_b="results/sigprofiler/de_novo/sigprofiler_W_de_novo.csv",
    output_path="results/analysis/cosine_similarity_matrix.csv",
)


matches = match_signatures(
    similarity_matrix=similarity_matrix,
    output_path="results/analysis/matched_signatures.csv",
)


# Save heatmap of pairwise cosine similarities
plot_similarity_heatmap(
    similarity_matrix=similarity_matrix,
    output_path="figures/cosine_similarity_heatmap.png",
)


# Example dot plot for matched signature similarities
matches_for_plot = matches.copy()
matches_for_plot["comparison"] = (
    matches_for_plot["signature_a"]
    + " vs "
    + matches_for_plot["signature_b"]
)
matches_for_plot["mutation_type"] = "SBS"
matches_for_plot["method"] = "MuSiCal vs SigProfiler"

plot_cosine_similarity_dotplot(
    similarity_df=matches_for_plot,
    output_path="figures/cosine_similarity_dotplot.png",
    title="Matched de novo signature similarity",
)


# Reconstruction quality
ev_df = compute_explained_variance_table(
    runs=[
        {
            "mutation_type": "SBS",
            "method": "MuSiCal",
            "observed_matrix": INPUT_MATRIX,
            "signature_matrix": "results/musical/de_novo/musical_W_de_novo.csv",
            "exposure_matrix": "results/musical/de_novo/musical_H_de_novo.csv",
        },
        {
            "mutation_type": "SBS",
            "method": "SigProfiler",
            "observed_matrix": INPUT_MATRIX,
            "signature_matrix": "results/sigprofiler/de_novo/sigprofiler_W_de_novo.csv",
            "exposure_matrix": "results/sigprofiler/de_novo/sigprofiler_H_de_novo.csv",
        },
    ],
    output_path="results/analysis/explained_variance.csv",
)


plot_explained_variance(
    ev_df=ev_df,
    output_path="figures/explained_variance.png",
)


print("Analysis completed.")
