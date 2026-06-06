import re
import pandas as pd
import matplotlib.pyplot as plt

# Input LME results file
path = "lme_results_all_group_plus_sex_fixed_with_mean_difference.csv"

# Output PDF
out_path = "p_value_histograms_by_analysis_mode_faceted.pdf"

# Read data
df = pd.read_csv(path)
df.columns = df.columns.str.strip()

# Required columns
metabolite_col = "metabolite"
p_col = "p_value_group"

required = {metabolite_col, p_col}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Coerce p-values to numeric
df[p_col] = pd.to_numeric(df[p_col], errors="coerce")

# Extract analysis mode from metabolite prefixes:
# C18_POS_, C18_NEG_, HILIC_POS_, HILIC_NEG_
def get_analysis_mode(metabolite):
    metabolite = str(metabolite)
    match = re.match(r"^(C18_POS|C18_NEG|HILIC_POS|HILIC_NEG)_", metabolite)
    if match:
        return match.group(1)
    return "UNKNOWN"

df["analysis_mode"] = df[metabolite_col].apply(get_analysis_mode)

# Drop invalid p-values
df = df.dropna(subset=[p_col]).copy()
df = df[(df[p_col] >= 0) & (df[p_col] <= 1)].copy()

# Desired facet order
mode_order = ["C18_POS", "C18_NEG", "HILIC_POS", "HILIC_NEG"]

# Prepare one figure with four facets
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=False)
fig.patch.set_facecolor("white")

axes = axes.flatten()

for ax, mode in zip(axes, mode_order):
    sub = df[df["analysis_mode"] == mode].copy()

    ax.set_facecolor("white")

    ax.hist(
        sub[p_col],
        bins=20,
        range=(0, 1),
        edgecolor="black"
    )

    ax.axvline(0.05, linestyle="--", linewidth=1)

    n_total = len(sub)
    n_sig = int((sub[p_col] < 0.05).sum())

    ax.set_title(f"{mode}\nn = {n_total}, p < 0.05: {n_sig}")
    ax.set_xlabel("p-value")
    ax.set_ylabel("Number of metabolites")

    ax.grid(False)

    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_linewidth(1.0)

# Remove any unused axes if fewer than 4 modes are present
for ax in axes[len(mode_order):]:
    ax.axis("off")

plt.tight_layout()

plt.savefig(out_path, format="pdf", bbox_inches="tight")
plt.close()

print(f"Wrote {out_path}")

# Print counts for sanity check
summary = (
    df.groupby("analysis_mode")
      .agg(
          n_metabolites=(p_col, "size"),
          n_p_lt_0_05=(p_col, lambda x: int((x < 0.05).sum()))
      )
      .reset_index()
)

print(summary.to_string(index=False))
