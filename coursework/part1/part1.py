import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize

def speedup(p,f):
    return 1/(f + (1-f)/p)

def efficiency(p,f):
    return speedup(p,f)/p


funcs = ["init","step","norm","dxdt",""]
variants = {}
variants["init"]= ["default","collapse","manual"] # no default
variants["step"]= ["default","collapse","manual"]
variants["norm"]= ["default","collapse","manual"] #no manual
variants["dxdt"]= ["default","collapse","single","combo"] # none yet
variants[""] = [""]
colors = ["red","green","blue","purple"]
for func in funcs:
    fig, axs = plt.subplots(1,2,figsize =(8,6))
    axs[0].set_xlabel("Number of cores")
    axs[0].set_ylabel("Speedup")
    axs[1].set_xlabel("Number of cores")
    axs[1].set_ylabel("Efficiency")
    if (func != ""):
        axs[0].set_title(f"Strong Scaling when parallelising {func}")
        axs[1].set_title(f"Efficiency when parallelising {func}")
    else:
        axs[0].set_title("Strong Scaling Speedup for overall paralellised program")
        axs[1].set_title("Strong scaling Efficiency for overall paralellised program")
    data  = {}
    for i,variant in enumerate(variants[func]):
        if (func != ""):
            filename = f"newData/part1_{func}_{variant}_strong_scaling.dat"
        else:
            filename = f"newData/part1_strong_scaling.dat"
        data[variant] = pd.read_csv(filename,delimiter=",")
        data[variant]["cores"] = data[variant]["cores"].astype(int)
        data[variant] = data[variant][data[variant]["cores"]<=64]
        data[variant]["time_seconds"] = pd.to_numeric(data[variant]["time_seconds"], errors="coerce")
        grouped = data[variant].groupby("cores")["time_seconds"]

        mean_data = grouped.agg(
            mean="mean",
            std="std",
            count="count"
        ).reset_index()
        mean_data["stderr"] = mean_data["std"] / np.sqrt(mean_data["count"])
        # Determine single-core time
        if 1 in mean_data["cores"].values:
            single_core_time = mean_data.loc[mean_data["cores"] == 1, "mean"].values[0]
        else:
            single_core_time = mean_data["mean"].iloc[0]  # fallback to smallest core count
        
        mean_data["speedup"] = single_core_time / mean_data["mean"]
        mean_data["speedup_error"] = (single_core_time / mean_data["mean"]**2) * mean_data["stderr"]
        mean_data["stderr"] = mean_data["std"] / np.sqrt(mean_data["count"])
        mean_data["efficiency"] = mean_data["speedup"] / mean_data["cores"]
        mean_data["efficiency_error"] = mean_data["speedup_error"] / mean_data["cores"]
        mean_data.to_csv(f"processedData/part1_{func}_{variant}_strong_scaling_processed.csv",index=False)
        xs = np.arange(1, max(mean_data["cores"])+1,0.1)
        params,covar = optimize.curve_fit(speedup,mean_data["cores"],mean_data["speedup"],p0=[0.1],sigma =mean_data["speedup_error"],absolute_sigma=True,bounds=(0,1))
        f = params[0]
        axs[0].plot(xs,speedup(xs,f), color=colors[i])
        axs[0].errorbar(mean_data["cores"],mean_data["speedup"],yerr=mean_data["speedup_error"],marker='',capsize=3,label = variant+f" (f={f:.2f})",color=colors[i],linestyle='')
        axs[1].plot(xs,efficiency(xs,f), color=colors[i])
        axs[1].errorbar(mean_data["cores"],mean_data["efficiency"],yerr=mean_data["efficiency_error"],marker='',capsize=3,label = variant+f" (f={f:.2f})",color=colors[i],linestyle='')
    axs[0].legend()
    axs[1].legend()
    fig.tight_layout()
    plt.savefig(f"plots/part1_{func}_strong_scaling.png")
    plt.show()
    