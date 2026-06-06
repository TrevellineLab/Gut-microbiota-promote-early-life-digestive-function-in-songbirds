import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the sex-adjusted LME results file with mean difference
path = "lme_results_all_group_plus_sex_fixed_with_mean_difference.csv"
df = pd.read_csv(path)

# Standardize column names
df.columns = df.columns.str.strip()

# Use unadjusted mean difference for x-axis and group p-value for y-axis
effect_col = "mean_difference_CV_minus_AX"
p_col = "p_value_group"

# Basic validation
required = {effect_col, p_col}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Coerce to numeric and remove invalid rows
df[effect_col] = pd.to_numeric(df[effect_col], errors="coerce")
df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
df = df.dropna(subset=[effect_col, p_col]).copy()

# Remove p = 0 or negative p-values because -log10 is undefined/infinite
df = df[df[p_col] > 0].copy()

# Compute -log10(p)
df["neg_log10_p"] = -np.log10(df[p_col])

# Determine significance
df["significant"] = df[p_col] < 0.05

# Prepare figure: single plot, white background, no gridlines, no title
fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Non-significant points
mask_ns = ~df["significant"]
ax.scatter(
    df.loc[mask_ns, effect_col],
    df.loc[mask_ns, "neg_log10_p"],
    s=18,
    alpha=0.7,
    label="p >= 0.05"
)

# Significant points with different marker shape
mask_sig = df["significant"]
ax.scatter(
    df.loc[mask_sig, effect_col],
    df.loc[mask_sig, "neg_log10_p"],
    s=24,
    alpha=0.85,
    label="p < 0.05"
)

# Reference lines
y_thresh = -np.log10(0.05)
ax.axhline(y_thresh, linestyle="--", linewidth=1)
ax.axvline(0, linewidth=1)

# Labels
ax.set_xlabel("Mean difference (CV - AX)")
ax.set_ylabel("-log10(p-value)")

# Axes styling
ax.grid(False)
for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_linewidth(1.0)

# Legend: bottom-right, boxed
leg = ax.legend(loc="lower right", frameon=True)
leg.get_frame().set_linewidth(1.0)

plt.tight_layout()

# Save as PDF
out_path = "fig.2.e.plot.pdf"
plt.savefig(out_path, format="pdf", bbox_inches="tight")
plt.close()

print(f"Wrote {out_path}")
print(f"Plotted {len(df)} metabolites; significant at p < 0.05: {int(df['significant'].sum())}")
