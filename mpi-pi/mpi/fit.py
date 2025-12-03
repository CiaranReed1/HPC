import matplotlib.pyplot as plt    
import numpy as np
import pandas as pd
import scipy.optimize as sp

def strong_model(p,f):
    """ Amdahl's law model for strong scaling
    p : number of processes
    f : fraction of serial code
    """
    return 1.0/ (f + (1-f)/p)

strong_data = pd.read_csv("pll_strong_ham.dat",delimiter=",")
strong_data["speedup"] = strong_data["time"].iloc[0]/strong_data["time"]
weights = 1 / strong_data["speedup"] 
param,covar = sp.curve_fit(strong_model,strong_data["size"],strong_data["speedup"],p0=[0.9],bounds=(0,1),sigma=weights)
fitted_f = param[0]
xs = np.arange(1,16,0.1)
ys = strong_model(xs,fitted_f)


fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(strong_data["size"],strong_data["speedup"],marker="x",linestyle = "")
ax.plot(xs, ys, label=f"Fitted model (f={fitted_f:.4f})",color="red")
#ax.set_xscale("log",base=2)
#ax.set_yscale("log",base=2)
ax.set_xlabel("Number of processes")
ax.set_ylabel("Speedup")
ax.set_title("Strong Scaling of MPI Pi Calculation (N = 10000000)")
ax.legend()
plt.show()
fig.savefig("strong_scaling_fit.png")