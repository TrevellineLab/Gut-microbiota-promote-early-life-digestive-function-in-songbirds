#This script uses the BH method to FDR correct the p-values derived from linear mixed effects models testing for differences between CV and AX birds. BH correction was performed within each annotation category and within each RPLC-MS/MS or HILIC-MS/MS analysis mode (HILIC_POS, HILIC_NEG, C18_POS, or C18_NEG). The uncorrected p-values are available in Table S7. 

import numpy as np
from statsmodels.stats.multitest import multipletests

# Your two sets
C18_NEG_BileAcids = np.array([
    0.006216703,
    0.008851555,
    0.008954321,
    0.01859532,
    0.093817422
])

HILIC_NEG_BileAcids = np.array([
    0.013372141,
    0.040149952,
    0.065778224
])

def bh_correct(pvals):
    # returns BH-adjusted p-values in original order
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
    return qvals

# Apply independently
q1 = bh_correct(set1)
q2 = bh_correct(set2)

print("Set 1 BH-corrected:")
print(q1)

print("\nSet 2 BH-corrected:")
print(q2)