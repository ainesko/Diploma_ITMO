# Import necessary modules
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import time
import scipy as sp
import pickle
import musical

#De novo signature discovery
#The input mutation count matrix X should be of size (n_features, n_samples), where n_features = 96 for SBS signatures.
X=pd.read_csv('/path/to/MuSiCal/examples/data/simulated_example.Skin.Melanoma.X.csv', index_col=0)
model = musical.DenovoSig(X, 
                          min_n_components=1, # Minimum number of signatures to test
                          max_n_components=15, # Maximum number of signatures to test
                          init='random', # Initialization method
                          method='mvnmf', # mvnmf or nmf
                          n_replicates=20, # Number of mvnmf/nmf replicates to run per n_components
                          ncpu=1, # Number of CPUs to use
                          max_iter=100000, # Maximum number of iterations for each mvnmf/nmf run
                          bootstrap=True, # Whether or not to bootstrap X for each run
                          tol=1e-8, # Tolerance for claiming convergence of mvnmf/nmf
                          verbose=1, # Verbosity of output
                          normalize_X=False # Whether or not to L1 normalize each sample in X before mvnmf/nmf
                         )
model.fit()

# Number of discovered de novo signatures
print(model.n_components)

#Save de novo signatures
W_1 = pd.DataFrame(model.W)
W_1.to_csv('W_SkinMelanoma.csv')

#Signature assignment
thresh_grid = np.array([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5.])
catalog = musical.load_catalog('COSMIC-MuSiCal_v3p2_SBS_WGS') # or 'COSMIC_v3p1_Indel' for indels
W_catalog = catalog.W
print(W_catalog.shape[1])
model.assign_grid(W_catalog, 
                  method_assign='likelihood_bidirectional', # Method for performing matching and refitting
                  thresh_match_grid=thresh_grid, # Grid of threshold for matchinng
                  thresh_refit_grid=thresh_grid, # Grid of threshold for refitting
                  thresh_new_sig=0.0, # De novo signatures with reconstructed cosine similarity below this threshold will be considered novel
                  connected_sigs=False, # Whether or not to force connected signatures to co-occur
                  clean_W_s=False # An optional intermediate step to avoid overfitting to small backgrounds in de novo signatures for 96-channel SBS signatures
                 )
  
#in silico validation
model.validate_grid(validate_n_replicates=1, # Number of simulation replicates to perform for each grid point
                    grid_selection_method='pvalue', # Method for selecting the best grid point
                    grid_selection_pvalue_thresh=0.05 # Threshold used for selecting the best grid point
                   )
#Results
print(model.best_grid_point)
print(model.thresh_match)
print(model.thresh_refit)
W_s = model.W_s
H_s = model.H_s
print(W_s.columns.tolist())
  
#Save final matched signatures and signature assignments
W_s  = pd.DataFrame(data=W_s[1:,1:], index=W_s [1:,0], columns=W_s [0,1:])
H_s  = pd.DataFrame(data=H_s[1:,1:], index=H_s [1:,0], columns=H_s [0,1:])
W_s.to_csv('W_s_SkinMelanoma.csv')
H_s.to_csv('H_s_SkinMelanoma.csv')
