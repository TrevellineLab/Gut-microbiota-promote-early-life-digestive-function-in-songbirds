import numpy as np
from scipy import stats

# Your paired data
x = np.array([
    0.433771943,
    0.282173258,
    0.217283618,
    0.237128772,
    0.210110815
])

y = np.array([
    -0.375935683,
    -0.244550155,
    -0.188312469,
    -0.205511601,
    -0.182096033
])

# Paired t-test
t_stat, p_value = stats.ttest_rel(x, y)

print("t-statistic:", t_stat)
print("p-value:", p_value)

##Results##
#Mean difference (x − y) = 0.5153748694
#SD of differences = 0.1726849
#n = 5, df = 4
#Paired t-test t-statistic = 6.6709; two-sided p-value = 0.00263