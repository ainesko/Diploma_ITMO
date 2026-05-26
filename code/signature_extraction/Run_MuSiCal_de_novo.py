"""Extract de novo mutational signatures using MuSiCal."""

from pathlib import Path
import pandas as pd
import musical


def run_musical_denovo(
    input_matrix: str | Path,
    output_dir: str | Path,
    sig_type: str = "SBS96",
    min_components: int = 1,
    max_components: int = 6,
    n_replicates: int = 20,
    ncpu: int = 1,
):
    """
    Run de novo mutational signature extraction with MuSiCal.

    Parameters
    ----------
    input_matrix : str or Path
        Path to a mutation count matrix.

        Expected format:
        - rows: mutation channels
        - columns: samples
        - values: mutation counts

        For SBS96, the matrix should contain 96 trinucleotide
        substitution contexts.

        For ID83, the matrix should contain 83 indel contexts.

    output_dir : str or Path
        Directory where the extracted signature and exposure
        matrices will be saved.

    sig_type : str
        Type of mutational context:
        - "SBS96" for single base substitutions
        - "ID83" for insertions and deletions

    min_components : int
        Minimum number of de novo signatures to test.

    max_components : int
        Maximum number of de novo signatures to test.

    n_replicates : int
        Number of repeated mvNMF runs for each tested number
        of signatures.

    ncpu : int
        Number of CPU cores used for computation.

    Returns
    -------
    W : pandas.DataFrame
        De novo signature matrix:
        mutation channels × extracted signatures.

    H : pandas.DataFrame
        Exposure matrix:
        samples × extracted signatures.

    model : musical.DenovoSig
        Fitted MuSiCal model.
    """

    input_matrix = Path(input_matrix)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load mutation count matrix.
    X = pd.read_csv(input_matrix, index_col=0)

    # Run de novo signature extraction.
    model = musical.DenovoSig(
        X,
        min_n_components=min_components,
        max_n_components=max_components,
        init="random",
        method="mvnmf",
        n_replicates=n_replicates,
        ncpu=ncpu,
        max_iter=100000,
        bootstrap=True,
        tol=1e-8,
        verbose=1,
        normalize_X=False,
    )

    model.fit()

    print(f"Number of discovered de novo signatures: {model.n_components}")

    # Save de novo signatures.
    W = pd.DataFrame(
        model.W,
        index=X.index,
        columns=[f"Signature_{i + 1}" for i in range(model.W.shape[1])],
    )

    # Save sample exposures.
    H = pd.DataFrame(
        model.H.T,
        index=X.columns,
        columns=[f"Signature_{i + 1}" for i in range(model.H.shape[0])],
    )

    W.to_csv(output_dir / "musical_W_de_novo.csv")
    H.to_csv(output_dir / "musical_H_de_novo.csv")

    # Optional quick visualization.
    if sig_type == "ID83":
        musical.sigplot_bar(W, sig_type="Indel83")
    else:
        musical.sigplot_bar(W)

    return W, H, model


if __name__ == "__main__":
    run_musical_denovo(
        input_matrix="data/example_matrix.csv",
        output_dir="results/musical/de_novo",
        sig_type="SBS96",
        min_components=1,
        max_components=6,
        n_replicates=20,
        ncpu=1,
    )
