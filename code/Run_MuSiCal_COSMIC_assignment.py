"""Assign previously extracted MuSiCal de novo signatures to COSMIC signatures."""

from pathlib import Path

import numpy as np
import pandas as pd
import musical


def run_musical_cosmic_assignment(
    signature_matrix: str | Path,
    output_dir: str | Path,
    sig_type: str = "SBS96",
):
    """
    Match previously extracted de novo signatures to COSMIC signatures.

    Parameters
    ----------
    signature_matrix : str or Path
        Path to a previously extracted de novo signature matrix.

        Expected format:
        - rows: mutation channels
        - columns: extracted signatures

        Example:
                            Signature_1   Signature_2
            A[C>A]A              0.01          0.03
            A[C>A]C              0.00          0.05

    output_dir : str or Path
        Directory where COSMIC assignment results will be saved.

    sig_type : str
        Type of mutational context:
        - "SBS96" for single base substitutions
        - "ID83" for insertions and deletions

    Returns
    -------
    W_s : pandas.DataFrame
        COSMIC-assigned signature matrix.
    """

    signature_matrix = Path(signature_matrix)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load previously extracted de novo signatures.
    W = pd.read_csv(signature_matrix, index_col=0)

    # Load appropriate COSMIC reference catalog.
    if sig_type == "ID83":
        catalog = musical.load_catalog("COSMIC_v3p1_Indel")
    else:
        catalog = musical.load_catalog("COSMIC-MuSiCal_v3p2_SBS_WGS")

    W_catalog = catalog.W

    # Threshold grid used for matching.
    thresh_grid = np.array([
        0.0001, 0.0002, 0.0005,
        0.001, 0.002, 0.005,
        0.01, 0.02, 0.05,
        0.1, 0.2, 0.5,
        1.0, 2.0, 5.0,
    ])

    # Match de novo signatures to COSMIC signatures.
    W_s, sig_map, cosine = musical.refit.match(
        W,
        W_catalog,
        method="likelihood_bidirectional",
        thresh=0.001,
    )

    print("Matched COSMIC signatures:")
    print(sig_map)

    print("\nCosine similarities:")
    print(cosine)

    # Save matched signatures.
    W_s.to_csv(output_dir / "musical_W_COSMIC_assigned.csv")

    return W_s


if __name__ == "__main__":
    run_musical_cosmic_assignment(
        signature_matrix="results/musical/de_novo/musical_W_de_novo.csv",
        output_dir="results/musical/cosmic_assignment",
        sig_type="SBS96",
    )
