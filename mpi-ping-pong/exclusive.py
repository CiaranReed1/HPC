import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  

data_onenode = pd.read_csv('exclusive-onenode.dat', delimiter=',')
data_twonode = pd.read_csv('exclusive-twonode.dat', delimiter=',')
data_fournode = pd.read_csv('exclusive-fournode.dat', delimiter=',')
fig, axs = plt.subplots(1, 3, figsize=(10, 8))
for i, data in enumerate([data_onenode, data_twonode,data_fournode]):
    axs[i].errorbar(data['nbytes'], data['mean'], yerr=data['stdev'],linestyle = "", label="nodes  = {}".format(i+1),color = "blue", capsize=5)
    axs[i].plot(data['nbytes'], data['min'], color = "grey",linestyle = "",marker = "X")
    axs[i].plot(data['nbytes'], data['max'], color = "grey",linestyle = "",marker = "X")
    axs[i].set_xlabel('nbytes')
    axs[i].set_xscale('log', base=2)
    axs[i].set_yscale('log', base=2)
    axs[i].set_ylabel('time (s)')
    axs[i].legend()
fig.savefig('exclusive.png')
plt.show()
print("success")