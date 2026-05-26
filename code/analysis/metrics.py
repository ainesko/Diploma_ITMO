"""Metrics for mutational signature analysis."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cosine


def cosine_similarity(x, y):
    """Compute cosine similarity between two vectors."""
    return 1 - cosine(x, y)


def compute_similarity_matrix(
    signatures_a,
    signatures_b,
    output_path=None,
):
    """
    Compute pairwise cosine similarity between two signature matrices.

    Parameters
    ----------
    signatures_a : str, Path, or pandas.DataFrame
        First signature matrix.

        Expected format:
        - rows: mutation channels
        - columns: signatures

    signatures_b : str, Path, or pandas.DataFrame
        Second signature matrix.

        Expected format:
        - rows: mutation channels
        - columns: signatures

    output_path : str or Path, optional
        Path where the cosine similarity matrix will be saved.

    Returns
    -------
    similarity_matrix : pandas.DataFrame
        Pairwise cosine similarity matrix.
    """

    if not isinstance(signatures_a, pd.DataFrame):
        signatures_a = pd.read_csv(signatures_a, index_col=0)

    if not isinstance(signatures_b, pd.DataFrame):
        signatures_b = pd.read_csv(signatures_b, index_col=0)

    common_channels = signatures_a.index.intersection(signatures_b.index)

    signatures_a = signatures_a.loc[common_channels]
    signatures_b = signatures_b.loc[common_channels]

    similarity_matrix = pd.DataFrame(
        index=signatures_a.columns,
        columns=signatures_b.columns,
        dtype=float,
    )

    for sig_a in signatures_a.columns:
        for sig_b in signatures_b.columns:
            similarity_matrix.loc[sig_a, sig_b] = cosine_similarity(
                signatures_a[sig_a].values,
                signatures_b[sig_b].values,
            )

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        similarity_matrix.to_csv(output_path)

    return similarity_matrix


def match_signatures(
    similarity_matrix,
    output_path=None,
):
    """
    Match signatures using the Hungarian algorithm.

    Parameters
    ----------
    similarity_matrix : str, Path, or pandas.DataFrame
        Pairwise cosine similarity matrix.

    output_path : str or Path, optional
        Path where the matched signature pairs will be saved.

    Returns
    -------
    matches : pandas.DataFrame
        Table with matched signatures and cosine similarities.
    """

    if not isinstance(similarity_matrix, pd.DataFrame):
        similarity_matrix = pd.read_csv(similarity_matrix, index_col=0)

    row_ind, col_ind = linear_sum_assignment(-similarity_matrix.values)

    matches = pd.DataFrame({
        "signature_a": similarity_matrix.index[row_ind],
        "signature_b": similarity_matrix.columns[col_ind],
        "cosine_similarity": similarity_matrix.values[row_ind, col_ind],
    })

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        matches.to_csv(output_path, index=False)

    return matches


def compute_explained_variance(
    observed_matrix,
    signature_matrix,
    exposure_matrix,
):
    """
    Compute explained variance for one reconstructed mutation matrix.

    Parameters
    ----------
    observed_matrix : str, Path, or pandas.DataFrame
        Original mutation count matrix.

        Expected format:
        - rows: mutation channels
        - columns: samples

    signature_matrix : str, Path, or pandas.DataFrame
        Signature matrix W.

        Expected format:
        - rows: mutation channels
        - columns: signatures

    exposure_matrix : str, Path, or pandas.DataFrame
        Exposure matrix H.

        Expected format:
        - rows: samples
        - columns: signatures

    Returns
    -------
    explained_variance : float
        Fraction of variance explained by the reconstructed matrix.
    """

    if not isinstance(observed_matrix, pd.DataFrame):
        observed_matrix = pd.read_csv(observed_matrix, index_col=0)

    if not isinstance(signature_matrix, pd.DataFrame):
        signature_matrix = pd.read_csv(signature_matrix, index_col=0)

    if not isinstance(exposure_matrix, pd.DataFrame):
        exposure_matrix = pd.read_csv(exposure_matrix, index_col=0)

    common_channels = observed_matrix.index.intersection(signature_matrix.index)
    common_samples = observed_matrix.columns.intersection(exposure_matrix.index)
    common_signatures = signature_matrix.columns.intersection(
        exposure_matrix.columns
    )

    X = observed_matrix.loc[common_channels, common_samples]
    W = signature_matrix.loc[common_channels, common_signatures]
    H = exposure_matrix.loc[common_samples, common_signatures]

    reconstructed = pd.DataFrame(
        np.dot(W.values, H.values.T),
        index=W.index,
        columns=H.index,
    )

    residual_sum_of_squares = np.sum(
        (X.values - reconstructed.values) ** 2
    )

    total_sum_of_squares = np.sum(
        (X.values - X.values.mean()) ** 2
    )

    explained_variance = 1 - residual_sum_of_squares / total_sum_of_squares

    return explained_variance


def compute_explained_variance_table(
    runs,
    output_path=None,
):
    """
    Compute explained variance for several datasets and methods.

    Parameters
    ----------
    runs : list of dict
        Each dictionary describes one reconstruction.

        Required keys:
        - mutation_type
        - method
        - observed_matrix
        - signature_matrix
        - exposure_matrix

        Example:
        {
            "mutation_type": "SBS",
            "method": "MuSiCal",
            "observed_matrix": "data/example_sbs_matrix.csv",
            "signature_matrix": "results/musical/de_novo/musical_W_de_novo.csv",
            "exposure_matrix": "results/musical/de_novo/musical_H_de_novo.csv",
        }

    output_path : str or Path, optional
        Path where the explained variance table will be saved.

    Returns
    -------
    ev_df : pandas.DataFrame
        Table with explained variance values.
    """

    records = []

    for run in runs:
        explained_variance = compute_explained_variance(
            observed_matrix=run["observed_matrix"],
            signature_matrix=run["signature_matrix"],
            exposure_matrix=run["exposure_matrix"],
        )

        records.append({
            "mutation_type": run["mutation_type"],
            "method": run["method"],
            "explained_variance": explained_variance,
        })

    ev_df = pd.DataFrame(records)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        ev_df.to_csv(output_path, index=False)

    return ev_df
