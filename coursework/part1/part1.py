import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize

def speedup(p,f):
    return 1/(f + (1-f)/p)


funcs = ["init","step","norm"]#,"dxdt"]
variants = {}
variants["init"]= ["collapse","manual"] # no default
variants["step"]= ["default","collapse","manual"]
variants["norm"]= ["default","collapse"] #no manual
#variants["dxdt"]= ["default","collapse","single"] # none yet
colors = ["red","green","blue"]
for func in funcs:
    fig, ax = plt.subplots(1,1,figsize =(8,6))
    ax.set_xlabel("Number of cores")
    ax.set_ylabel("Speedup")
    ax.set_title(f"Strong Scaling for {func} function")
    data  = {}
    for i,variant in enumerate(variants[func]):
        filename = f"newData/part1_{func}_{variant}_strong_scaling.dat"
        data[variant] = pd.read_csv(filename,delimiter=",")
        data[variant]["cores"] = data[variant]["cores"].astype(int)
        data[variant]["time_seconds"] = pd.to_numeric(data[variant]["time_seconds"], errors="coerce")
        mean_data = data[variant].groupby("cores", as_index=False)["time_seconds"].mean()
        mean_data["stderr"] = data[variant].groupby("cores")["time_seconds"].std()/np.sqrt(data[variant].groupby("cores")["time_seconds"].count()).values
        
        # Determine single-core time
        if 1 in mean_data["cores"].values:
            single_core_time = mean_data.loc[mean_data["cores"] == 1, "time_seconds"].values[0]
        else:
            single_core_time = mean_data["time_seconds"].iloc[0]  # fallback to smallest core count
        
        mean_data["speedup"] = single_core_time / mean_data["time_seconds"]
        mean_data["speedup_error"] = (single_core_time / mean_data["time_seconds"]**2) * mean_data["stderr"]
        
        xs = np.arange(1, max(mean_data["cores"])+1,0.1)
        params,covar = optimize.curve_fit(speedup,mean_data["cores"],mean_data["speedup"],p0=[0.1])
        f = params[0]
        ax.plot(xs,speedup(xs,f), color=colors[i])
        ax.errorbar(mean_data["cores"],mean_data["speedup"],yerr=mean_data["speedup_error"],marker='x',label = variant+f" (f={f:.2f})",color=colors[i],linestyle='')
    ax.legend()
    plt.savefig(f"plots/part1_{func}_strong_scaling.png")
    plt.show()
    