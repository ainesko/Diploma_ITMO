import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity

# Load the signatures extracted by MuSiCal method
MSC_After_B_Post = pd.read_csv('/Users/ksenia/Desktop/ITMO/MSC_W_After_B_Post_100.csv', index_col=0)
# Load the signatures extracted by SigProfiler method
SP_After_B_Post = pd.read_csv('/Users/ksenia/Desktop/ITMO/After_B_post_results_SigProfiler/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt', sep='\t', index_col=0)

# Reindex both DataFrames to follow the standard SBS96 COSMIC mutation type order
# (SBS96_COSMIC_ORDER should be defined elsewhere in your code)
MSC_After_B_Post = MSC_After_B_Post.reindex(SBS96_COSMIC_ORDER)
SP_After_B_Post = SP_After_B_Post.reindex(SBS96_COSMIC_ORDER)

# Compute cosine similarity between all pairs of signatures from the two methods
cos_sim = cosine_similarity(MSC_After_B_Post.T, SP_After_B_Post.T)

# Get number of signatures from each method
n_music = MSC_After_B_Post.shape[1]
n_sigprof = SP_After_B_Post.shape[1]

# Create descriptive labels for the heatmap axes
cols_matrix1 = [f"MSC_Sig_{i+1}" for i in range(n_music)]
cols_matrix2 = [f"SP_Sig_{i+1}" for i in range(n_sigprof)]

# Create figure for visualization
plt.figure(figsize=(10, 8))

# Generate heatmap using seaborn
sns.heatmap(
    cos_sim,                    # Cosine similarity matrix
    annot=True,                 # Display values in cells
    annot_kws={"size": 30, "weight": "bold", "color": "white"},  # Format for annotations
    fmt='.2f',                  # Format: 2 decimal places
    cmap="RdBu_r",              # Diverging colormap (red-blue reversed)
    center=0,                   # Center color scale at 0
    vmin=-1, vmax=1,            # Full range for cosine similarity
    xticklabels=cols_matrix2,   # SigProfiler signatures on x-axis
    yticklabels=cols_matrix1,   # MuSiCal signatures on y-axis
    square=True,                # Make cells square for better proportions
    cbar_kws={
        'label': 'Cosine Similarity', 
        'shrink': 0.75,
        'extend': 'both'        # Add arrows to colorbar ends
    },
    linewidths=0.5,             # Thin grid lines between cells
    linecolor='gray'            # Grid line color
)

# Set title with clear description of comparison
plt.title("MuSiCal vs SigProfiler De Novo Signatures\nCosine Similarity, After_B_Post (SBS)", 
          fontsize=20, pad=25, weight='bold')
plt.xlabel("SigProfiler Signatures", fontsize=18, labelpad=20)
plt.ylabel("MuSiCal Signatures", fontsize=18, labelpad=20)

# Format axis tick labels
plt.xticks(fontsize=16, rotation=0)   # Horizontal x-axis labels
plt.yticks(fontsize=16, rotation=90)  # Vertical y-axis labels

# Format colorbar
cbar = plt.gcf().axes[-1]
cbar.tick_params(labelsize=16)
cbar.set_ylabel('Cosine Similarity', fontsize=16, labelpad=15)

# Turn off grid (not needed for heatmap)
plt.grid(False)

# Adjust layout to prevent clipping
plt.tight_layout()

# Save figure as SVG vector graphic for publication quality
plt.savefig('/Users/ksenia/Desktop/ITMO/heatmap_SBS_B_Post.svg', 
            format='svg', 
            dpi=300,                     # High resolution
            bbox_inches='tight',         # Trim white space
            facecolor='white',           # Background color
            edgecolor='none')            # No border

# Display the plot
plt.show()
