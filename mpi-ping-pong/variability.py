import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  

data_onenode = pd.read_csv('variability-onenode.dat', delimiter=',')
data_twonode = pd.read_csv('variability-twonode.dat', delimiter=',')
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
colors = ["blue","red"]
for i, data in enumerate([data_onenode, data_twonode]):
    ax.errorbar(data['nbytes'], data['mean'], yerr=data['stdev'],linestyle = "", label="nodes  = {}".format(i+1),color = colors[i], capsize=5)
    ax.plot(data['nbytes'], data['min'], color = colors[i],linestyle = "",marker = "X")
    ax.plot(data['nbytes'], data['max'], color = colors[i],linestyle = "",marker = "X")
    ax.set_xlabel('nbytes')
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    ax.set_ylabel('time (s)')
    ax.legend()
fig.savefig('variability.png')
plt.show()
print("success")