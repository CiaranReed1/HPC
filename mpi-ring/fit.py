import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

data= pd.read_csv('ring.dat',delimiter=",")
fig, ax = plt.subplots(1,1,figsize =(8,6))
ax.plot(data["size"], data["mean_time"], marker='x',linestyle = "",color = "blue",label = "measured")
lr = linregress(data["size"], data["mean_time"])
x = np.arange(data["size"].min(), data["size"].max(),0.1)
ax.plot(x, lr.slope*x + lr.intercept, color = "red", label = f"fit: y={lr.slope:.2e}x + {lr.intercept:.2e}")
ax.set_xlabel("Number of processes")
ax.set_ylabel("Mean Runtime (s)")
ax.set_title("Mean Runtime of MPI Ring Reduction vs Number of Processes")
ax.legend()
fig.savefig("mpi_ring_performance.png", dpi=300)
plt.show()
