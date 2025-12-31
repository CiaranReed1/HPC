import numpy as np
import pandas as pd
from scipy import optimize
import matplotlib.pyplot as plt

def speedup_model(p,f):
    return f + (1-f)*p

data = pd.read_csv('scalingData/mpi-3_weak_scaling.dat',delimiter=",")
#print(data.head())
group_cols = ['problem_size', 'size']
initial_grouped = (
    data
    .groupby(group_cols)
    .agg(
        count=('init_time', 'size'),   # number of rows per group
        init_time_mean=('init_time', 'mean'),
        init_time_std=('init_time', 'std'),
        step_time_mean=('step_time', 'mean'),
        step_time_std=('step_time', 'std'),
        dxdt_time_mean=('dxdt_time', 'mean'),
        dxdt_time_std=('dxdt_time', 'std'),
        norm_time_mean=('norm_time', 'mean'),
        norm_time_std=('norm_time', 'std'),
        total_time_mean=('total_time', 'mean'),
        total_time_std=('total_time', 'std'),
    )
    .reset_index()
)
initial_grouped['init_time_se'] = initial_grouped['init_time_std'] / np.sqrt(initial_grouped['count'])
initial_grouped['step_time_se'] = initial_grouped['step_time_std'] / np.sqrt(initial_grouped['count'])
initial_grouped['dxdt_time_se'] = initial_grouped['dxdt_time_std'] / np.sqrt(initial_grouped['count'])
initial_grouped['norm_time_se'] = initial_grouped['norm_time_std'] / np.sqrt(initial_grouped['count'])
initial_grouped['total_time_se'] = initial_grouped['total_time_std'] / np.sqrt(initial_grouped['count'])
print(initial_grouped)

# Columns for which we want speedup
time_cols = ['init_time', 'step_time', 'dxdt_time', 'norm_time', 'total_time']

# Extract baseline (size=1) and max-size rows
baseline = initial_grouped[initial_grouped['size'] == 1].set_index('problem_size')
max_size_rows = initial_grouped.loc[initial_grouped.groupby('problem_size')['size'].idxmax()].set_index('problem_size')

# Initialize speedup dataframe
speedup_df = pd.DataFrame(index=baseline.index)

for col in time_cols:
    # Speedup
    speedup_df[f'{col}_speedup'] = baseline[f'{col}_mean'] / max_size_rows[f'{col}_mean']
    
    # Standard error propagation
    S = speedup_df[f'{col}_speedup']
    sigma1 = baseline[f'{col}_se']
    sigman = max_size_rows[f'{col}_se']
    speedup_df[f'{col}_speedup_se'] = S * np.sqrt((sigma1 / baseline[f'{col}_mean'])**2 + (sigman / max_size_rows[f'{col}_mean'])**2)

# Add max_size column
speedup_df['max_size'] = max_size_rows['size']

# Reset index to make problem_size a column
speedup_df = speedup_df.reset_index()

# Now safely set SE=0 for rows where max_size==1
se_cols = [f'{col}_speedup_se' for col in time_cols]
speedup_df.loc[speedup_df['max_size'] == 1, se_cols] = 0.0

print(speedup_df)
for col in time_cols:
    xs = np.arange(1, speedup_df['max_size'].max() + 1)
    popt, pcov = optimize.curve_fit(speedup_model, speedup_df['max_size'], speedup_df[f'{col}_speedup'], sigma=speedup_df[f'{col}_speedup_se'], absolute_sigma=True,p0=[0.9])
    ys = speedup_model(xs, *popt)
    yerr = np.sqrt(np.diag(pcov))
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    ax.errorbar(
        speedup_df['problem_size'],
        speedup_df[f'{col}_speedup'],
        yerr=speedup_df[f'{col}_speedup_se'],
        capsize=5,linestyle="",marker="",label="Measured Speedup",
        color = "blue"
    )
    ax.plot(xs, ys, label=f"Fitted Model: f={popt[0]:.3f} ± {yerr[0]:.3f}", color="red")
    ax.set_xlabel("Problem Size, n")
    ax.set_ylabel("Speedup")
    ax.set_title(f"Speedup for {col.replace('_', ' ').title()}")
    ax.legend()
    fig.savefig(f"plots/{col}_speedup.png",dpi=300)
