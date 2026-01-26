# Analysis of Somatic Mutational Signatures in Pediatric Cancers
## Authors
K.Mardanova¹, M.Artomov² 
1. ITMO University
2. Institute for Genomic Medicine Nationwide Children’s Hospital
## Project description

This study benchmarks mutational signature extraction methods on pediatric neuroblastoma data and evaluates adult-derived reference signatures.


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


## References

Alexandrov, L. B. et al. The repertoire of mutational signatures in human cancer. Nature 578, 94–101 (2020).
https://github.com/AlexandrovLab/SigProfilerExtractor

Jin, H. et al. Accurate and sensitive mutational signature analysis with MuSiCal. Nature Genetics, 56(3), 541–552 (2024). 
https://github.com/parklab/MuSiCal
