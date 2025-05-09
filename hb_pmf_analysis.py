# 0. Install dependencies if not running in a Google Colab Environment
# !pip install pandas numpy matplotlib scikit-learn scipy

# 1. Imports
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy.interpolate import interp1d
from math import comb

# 2. Auto‐detect CSV files
search_dirs = ['.', '/mnt/data']
csv_files = []
for base in search_dirs:
    if os.path.isdir(base):
        for root, _, files in os.walk(base):
            for f in files:
                if f.lower().endswith('.csv'):
                    csv_files.append(os.path.join(root, f))

# Identify required CSVs
sat_file   = next(f for f in csv_files if 'saturation' in os.path.basename(f).lower())
k_file     = next(f for f in csv_files if 'k values'   in os.path.basename(f).lower())
blood_file = next(f for f in csv_files if 'blood sample' in os.path.basename(f).lower())

# 3. Load data
blood_df = pd.read_csv(blood_file, encoding='utf-8')
k_df     = pd.read_csv(k_file, encoding='utf-8')
sat_df   = pd.read_csv(sat_file, encoding='utf-8')

# 4. Extract P50 & Hill coefficient h (row‐wise)
bp = blood_df.set_index(blood_df.columns[0])
p50_rows = [r for r in bp.index if 'p50' in r.lower() or 'p₅₀' in r.lower()]
h_rows   = [r for r in bp.index if 'nmax' in r.lower()]
P50 = bp.loc[p50_rows[0]].astype(float).mean()
h   = bp.loc[h_rows[0]].astype(float).mean()

# 5. Prepare saturation data
sat_df = sat_df.rename(columns={'pO2 (mmHg)':'pO2','Saturation':'S_raw'})
sat_df = sat_df[['pO2','S_raw']]
pO2_full = sat_df['pO2'].values
S_full   = sat_df['S_raw'].values

# 6. Parse K‐sets and keep only Winslow
key_col = 'Dataset' if 'Dataset' in k_df.columns else 'Paper'
all_sets = {
    row[key_col]: [
        float(str(row['K1']).split('+')[0]),
        float(str(row['K2']).split('+')[0]),
        float(str(row['K3']).split('+')[0]),
        float(str(row['K4']).split('+')[0]),
    ]
    for _, row in k_df.iterrows()
}
# Retain only the 'Winslow' entry
K_sets = {'Winslow': all_sets['Winslow']}

# 7. Ligand conversion
alpha = 0.0031    # mM O₂ per mmHg at 37 °C
use_alpha = False

# 8. PMF & helper functions
def pmf_binomial(pO2):
    p = pO2**h / (P50**h + pO2**h)
    return [comb(4, k) * p**k * (1-p)**(4-k) for k in range(5)]

def pmf_adair(ligand, K):
    cum = [1.0]
    for Ki in K:
        cum.append(cum[-1] * Ki)
    wts = [cum[i] * ligand**i for i in range(5)]
    tot = sum(wts)
    return [w/tot for w in wts]

def mean_sat(pmf):
    return sum(i * pmf[i] for i in range(5)) / 4

# 9. Compute curves & save results
results = {'pO2': pO2_full, 'S_raw': S_full}
results['Binomial'] = [mean_sat(pmf_binomial(p)) for p in pO2_full]
results['Winslow']  = [mean_sat(pmf_adair((alpha*p if use_alpha else p), K_sets['Winslow'])) for p in pO2_full]

res_df = pd.DataFrame(results)
res_df.to_csv("pmf_results.csv", index=False)

# 10. RMSE comparison
print("RMSE vs. Winslow data:")
for model in ['Binomial', 'Winslow']:
    rmse = np.sqrt(mean_squared_error(res_df['S_raw'], res_df[model]))
    print(f"{model:<10s} RMSE = {rmse:.4f}")

# 11. Validation at specific pO₂
levels = [10, 40, 60, 100]
interp = interp1d(pO2_full, S_full, fill_value="extrapolate")
validation = {'pO2': levels, 'Binomial': [], 'Binomial_Dev%': [], 'Winslow': [], 'Winslow_Dev%': []}

for p in levels:
    true_S = float(interp(p))
    bS = mean_sat(pmf_binomial(p))
    aS = mean_sat(pmf_adair((alpha*p if use_alpha else p), K_sets['Winslow']))
    validation['Binomial'].append(bS)
    validation['Binomial_Dev%'].append((bS-true_S)/true_S*100)
    validation['Winslow'].append(aS)
    validation['Winslow_Dev%'].append((aS-true_S)/true_S*100)

val_df = pd.DataFrame(validation)
val_df.to_csv("validation_specific_pO2.csv", index=False)

# 12. Plot mean‐saturation curves
plt.figure(figsize=(8,6))
plt.scatter(pO2_full, S_full, s=20, alpha=0.6, label='Winslow data')
plt.plot(res_df['pO2'], res_df['Binomial'], '--', label='Binomial')
plt.plot(res_df['pO2'], res_df['Winslow'],  label='Adair (Winslow)')
plt.xscale('log'); plt.xlabel('pO₂ (mmHg)'); plt.ylabel('Saturation')
plt.legend(); plt.tight_layout()
plt.savefig("pmf_comparison.pdf"); plt.show()

# 13. Bar‐plots at selected pO₂
fig, axs = plt.subplots(2, 2, figsize=(8, 6))
axs = axs.ravel()
for ax, p in zip(axs, levels):
    # Binomial
    pmfb = pmf_binomial(p)
    ax.bar(np.arange(5) - 0.3, pmfb, 0.2, label='Binomial')
    # Winslow
    pmfa = pmf_adair((alpha * p if use_alpha else p), K_sets['Winslow'])
    ax.bar(np.arange(5) + 0.2, pmfa, 0.2, label='Winslow')
    ax.set_xticks(range(5))
    ax.set_title(f"pO₂ = {p} mmHg")
    ax.set_xlabel("Bound O₂ (k)")
    ax.set_ylabel("P(k)")

# Collect handles/labels
handles, labels = axs[0].get_legend_handles_labels()
# Place legend at the bottom center
fig.legend(handles, labels,
           loc='lower center',
           bbox_to_anchor=(0.5, -0.05),
           ncol=2,
           fontsize=8)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("pmf_distributions.pdf")
print("✓ Saved pmf_distributions.pdf")
plt.show()

# 14. Line‐plots of P(k) vs pO₂ with consistent colors per k
plt.figure(figsize=(8,6))

# Use the first 5 colors from the default cycle
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color'][:5]

for k in range(5):
    c = color_cycle[k]
    # Binomial (dashed)
    probs_b = [pmf_binomial(p)[k] for p in pO2_full]
    plt.plot(pO2_full, probs_b, linestyle='--', color=c,
             label=f'Binomial, k={k}')
    # Winslow (solid)
    probs_w = [pmf_adair((alpha*p if use_alpha else p), K_sets['Winslow'])[k]
               for p in pO2_full]
    plt.plot(pO2_full, probs_w, linestyle='-', color=c,
             label=f'Winslow, k={k}')

plt.xscale('log')
plt.xlabel('pO₂ (mmHg)')
plt.ylabel('P(k)')
plt.title('PMF Profiles for Each k')
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left', fontsize=6)
plt.tight_layout()
plt.savefig("pmf_curves_each_state_matched_colors.pdf")
plt.show()

# 15. Sensitivity analysis on Winslow K-values
errors = [10, 18, 58, 52]
K_w = K_sets['Winslow']
sens = {
    'Winslow_low':  [mean_sat(pmf_adair((alpha*p if use_alpha else p), [K*(1-err/100) for K,err in zip(K_w,errors)])) for p in pO2_full],
    'Winslow_high': [mean_sat(pmf_adair((alpha*p if use_alpha else p), [K*(1+err/100) for K,err in zip(K_w,errors)])) for p in pO2_full]
}
sens_df = pd.concat([res_df, pd.DataFrame(sens)], axis=1)
sens_df.to_csv("pmf_results_sensitivity.csv", index=False)
print("✓ Written pmf_results_sensitivity.csv")

# Compute point-wise % deviation at every data point
residuals = 100*(res_df['Winslow'] - res_df['S_raw'])/res_df['S_raw']

# Plot residuals vs pO₂
plt.figure(figsize=(7,4))
plt.axhline(0, color='gray', linewidth=1, linestyle='--')
plt.scatter(res_df['pO2'], residuals, s=20, alpha=0.7)
plt.xscale('log')
plt.xlabel('pO₂ (mmHg)')
plt.ylabel('Deviation (%)')
plt.title('Residuals: Adair PMF (Winslow) vs. Actual Saturation')
plt.tight_layout()
plt.savefig("pmf_residuals_vs_pO2.pdf")
plt.show()
