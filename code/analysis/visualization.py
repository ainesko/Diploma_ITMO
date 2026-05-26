"""Visualization functions for mutational signature analysis."""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns


SBS_COLOR = "#163fa3"
ID_COLOR = "#7b2cbf"

LIGHT_SBS = "#dfe8ff"
LIGHT_ID = "#f0e2ff"

DARK_BLUE = "#0b2265"


def plot_similarity_heatmap(
    similarity_matrix,
    output_path,
    title="Cosine similarity between de novo signatures",
):
    """Plot heatmap of pairwise cosine similarities."""

    if not isinstance(similarity_matrix, pd.DataFrame):
        similarity_matrix = pd.read_csv(similarity_matrix, index_col=0)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(5.6, 4.8))

    sns.heatmap(
        similarity_matrix,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        vmin=0,
        vmax=1,
        square=True,
        linewidths=0.5,
        cbar_kws={"label": "Cosine similarity"},
    )

    plt.title(title, fontsize=13, fontweight="bold")
    plt.xlabel("Signature set B")
    plt.ylabel("Signature set A")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_explained_variance(
    ev_df,
    output_path,
    title="Reconstruction quality by mutation type",
):
    """
    Plot explained variance by mutation type and method.

    Expected columns:
    - mutation_type
    - method
    - explained_variance
    """

    if not isinstance(ev_df, pd.DataFrame):
        ev_df = pd.read_csv(ev_df)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ev_summary = (
        ev_df
        .groupby(["mutation_type", "method"])
        .agg(
            mean=("explained_variance", "mean"),
            std=("explained_variance", "std"),
            n=("explained_variance", "count"),
        )
        .reset_index()
    )

    ev_summary["sem"] = ev_summary["std"] / np.sqrt(ev_summary["n"])
    ev_summary["sem"] = ev_summary["sem"].fillna(0)

    mutation_order = ["SBS", "ID"]
    method_order = ["MuSiCal", "SigProfiler"]

    colors = {
        "SBS": {
            "edge": SBS_COLOR,
            "fill": LIGHT_SBS,
        },
        "ID": {
            "edge": ID_COLOR,
            "fill": LIGHT_ID,
        },
    }

    x = np.arange(len(mutation_order))
    width = 0.26

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for xi in x:
        panel = patches.FancyBboxPatch(
            (xi - 0.42, 0.0),
            0.84,
            1.03,
            boxstyle="round,pad=0.02,rounding_size=0.04",
            linewidth=0,
            facecolor="#f6f8ff",
            alpha=0.85,
            transform=ax.transData,
            zorder=0,
            clip_on=False,
        )
        ax.add_patch(panel)

    for i, method in enumerate(method_order):
        offset = (i - 0.5) * width

        means = []
        sems = []
        fill_colors = []
        edge_colors = []

        for mutation_type in mutation_order:
            row = ev_summary[
                (ev_summary["mutation_type"] == mutation_type)
                & (ev_summary["method"] == method)
            ].iloc[0]

            means.append(row["mean"])
            sems.append(row["sem"])
            fill_colors.append(colors[mutation_type]["fill"])
            edge_colors.append(colors[mutation_type]["edge"])

        bars = ax.bar(
            x + offset,
            means,
            width=width,
            yerr=sems,
            capsize=4,
            label=method,
            color=fill_colors,
            edgecolor=edge_colors,
            linewidth=2,
            hatch="///" if method == "SigProfiler" else None,
            alpha=0.6 if method == "SigProfiler" else 1,
            error_kw={
                "elinewidth": 1.2,
                "capthick": 1.2,
                "ecolor": "#202020",
            },
            zorder=3,
        )

        for bar, mean_val, sem_val in zip(bars, means, sems):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                mean_val + sem_val + 0.02,
                f"{mean_val:.2f}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                color="#1f1f1f",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        mutation_order,
        fontsize=12,
        fontweight="bold",
        color="#1f1f1f",
    )

    ax.set_ylabel(
        "Explained variance (mean ± SEM)",
        fontsize=11,
        color="#1f1f1f",
    )

    ax.set_ylim(0, 1.08)

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
        color=DARK_BLUE,
        pad=12,
    )

    legend = ax.legend(
        frameon=True,
        fontsize=10,
        loc="upper right",
    )

    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#d7dcf0")
    legend.get_frame().set_linewidth(1.0)

    ax.grid(axis="y", linestyle="--", linewidth=0.8, alpha=0.3)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_cosine_similarity_dotplot(
    similarity_df,
    output_path,
    x_col="comparison",
    y_col="cosine_similarity",
    group_col="mutation_type",
    method_col="method",
    title="Cosine similarity between extracted signatures",
):
    """
    Plot cosine similarity values as a grouped dot plot.

    This function is an example of the dot plots used for signature
    comparison. The same structure can be adapted for other comparison
    types, such as MuSiCal vs SigProfiler, pre- vs post-treatment,
    or arm-specific vs full-cohort signatures.

    Expected columns:
    - comparison
    - cosine_similarity
    - mutation_type
    - method
    """

    if not isinstance(similarity_df, pd.DataFrame):
        similarity_df = pd.read_csv(similarity_df)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mutation_order = ["SBS", "ID"]
    methods = list(similarity_df[method_col].dropna().unique())
    comparisons = list(similarity_df[x_col].dropna().unique())

    colors = {
        "SBS": SBS_COLOR,
        "ID": ID_COLOR,
    }

    markers = {
        method: marker
        for method, marker in zip(methods, ["o", "s", "D", "^"])
    }

    x = np.arange(len(comparisons))

    fig, ax = plt.subplots(figsize=(max(6, len(comparisons) * 0.8), 4.4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for i, method in enumerate(methods):
        method_data = similarity_df[
            similarity_df[method_col] == method
        ]

        offset = (i - (len(methods) - 1) / 2) * 0.12

        for mutation_type in mutation_order:
            subset = method_data[
                method_data[group_col] == mutation_type
            ]

            if subset.empty:
                continue

            x_positions = [
                comparisons.index(label) + offset
                for label in subset[x_col]
            ]

            ax.scatter(
                x_positions,
                subset[y_col],
                s=80,
                marker=markers.get(method, "o"),
                color=colors.get(mutation_type, "#4a4a4a"),
                edgecolor="#1f1f1f",
                linewidth=0.7,
                alpha=0.9,
                label=f"{mutation_type} {method}",
                zorder=3,
            )

            for xpos, value in zip(x_positions, subset[y_col]):
                ax.text(
                    xpos,
                    value + 0.025,
                    f"{value:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color="#1f1f1f",
                )

    ax.set_xticks(x)
    ax.set_xticklabels(
        comparisons,
        rotation=35,
        ha="right",
        fontsize=10,
    )

    ax.set_ylabel("Cosine similarity", fontsize=11)
    ax.set_ylim(0, 1.05)

    ax.set_title(
        title,
        fontsize=14,
        fontweight="bold",
        color=DARK_BLUE,
        pad=12,
    )

    ax.grid(axis="y", linestyle="--", linewidth=0.8, alpha=0.3)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))

    ax.legend(
        unique.values(),
        unique.keys(),
        frameon=True,
        fontsize=9.5,
        loc="lower right",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
