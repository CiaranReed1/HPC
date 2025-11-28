import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  

data_onenode = pd.read_csv('exclusive-onenode.dat', delimiter=',')
data_twonode = pd.read_csv('exclusive-twonode.dat', delimiter=',')
data_fournode = pd.read_csv('exclusive-fournode.dat', delimiter=',')
fig, axs = plt.subplots(1, 3, figsize=(10, 8))
fig2,axs2 = plt.subplots(1, 1, figsize=(10, 8))
colors = ['blue', 'orange', 'green']
for i, data in enumerate([data_onenode, data_twonode,data_fournode]):
    axs[i].errorbar(data['nbytes'], data['mean'], yerr=data['stdev'],linestyle = "", label="nodes  = {}".format(2**i),color = "blue", capsize=5)
    axs[i].plot(data['nbytes'], data['min'], color = "grey",linestyle = "",marker = "X")
    axs[i].plot(data['nbytes'], data['max'], color = "yellow",linestyle = "",marker = "X")
    axs[i].set_xlabel('nbytes')
    axs[i].set_xscale('log', base=2)
    axs[i].set_yscale('log', base=2)
    axs[i].set_ylabel('time (s)')
    axs[i].legend()
    axs2.set_xlabel('nbytes')
    axs2.set_xscale('log', base=2)
    axs2.set_yscale('log', base=2)
    axs2.set_ylabel('time (s)')
    axs2.plot(data['nbytes'], data['min'], color = colors[i],linestyle = "",marker = "X")
    axs2.plot(data['nbytes'], data['max'], color = colors[i],linestyle = "",marker = "d")
    axs2.errorbar(data['nbytes'], data['mean'], yerr=data['stdev'],linestyle = "", label="nodes  = {}".format(2**i), capsize=5,color = colors[i])
axs2.legend()
fig.savefig('exclusive.png')
fig2.savefig('exclusive-combined.png')
plt.show()
print("success")