# Analysis of Somatic Mutational Signatures in Pediatric Cancers
## Authors
K.Mardanova¹, M.Artomov² 
1. ITMO University
2. Institute for Genomic Medicine, Nationwide Children’s Hospital


## Project description

This project investigates somatic mutational signatures in pediatric neuroblastoma using de novo mutational signature extraction methods.

The analysis compares MuSiCal and SigProfilerExtractor, evaluates mutational profile stability between pre- and post-treatment samples, and explores limitations of adult-derived COSMIC reference signatures in pediatric cancers with low mutational burden.

## Installation of Python dependencies:

```bash
pip install -r requirements.txt
```

## MuSiCal installation

```bash
# Clone github repository
git clone https://github.com/parklab/MuSiCal

# Create a clean conda environment (recommended)
source ~/miniconda3/bin/activate root
conda create -n python37_musical python=3.7
conda install numpy scipy scikit-learn matplotlib pandas seaborn

# Install MuSiCal and dependencies
cd  /Path/To/MuSiCal
pip install ./MuSiCal
```

## SigProfilerExtractor installation
```bash
# Create a clean conda environment (recommended)
conda create -n sigprofiler python=3.9
conda activate sigprofiler

# Install SigProfilerExtractor and dependencies
pip install SigProfilerExtractor

# Install Reference Genome
from SigProfilerMatrixGenerator import install as gen_install
gen_install.install('GRCh38') 
```


## Workflow

1. De novo signature extraction with MuSiCal  
2. De novo signature extraction with SigProfilerExtractor  
3. COSMIC-based signature interpretation  
4. Comparison of extracted signatures  
5. Visualization and downstream analysis  

## Running the analysis

```bash
python code/run_analysis.py
```

## Input data

The pipeline expects mutation count matrices as input. Example format of an SBS96 mutation count matrix used as input for mutational signature extraction.

<img width="468" height="269" alt="image" src="https://github.com/user-attachments/assets/7941a6e4-4d2d-4021-a240-4e74d338c78c" />

Supported mutation contexts:

- SBS96
- ID83

For MuSiCal, input matrices should be provided in CSV format.

For SigProfilerExtractor, input matrices should be provided in TSV format.


## References

Alexandrov, L. B. et al. The repertoire of mutational signatures in human cancer. Nature 578, 94–101 (2020).
https://github.com/AlexandrovLab/SigProfilerExtractor

Jin, H. et al. Accurate and sensitive mutational signature analysis with MuSiCal. Nature Genetics, 56(3), 541–552 (2024). 
https://github.com/parklab/MuSiCal
