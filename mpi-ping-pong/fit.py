import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
data_1_node = pd.read_csv('one-node-pp.dat', delimiter=",", header=None, names=['iterations', 'nbytes','time'])
data_2_node = pd.read_csv('two-node-pp.dat', delimiter=",", header=None, names=['iterations', 'nbytes','time'])
fig, axs = plt.subplots(1, 2, figsize=(15, 6))
for i, data in enumerate([data_1_node, data_2_node]):
    axs[i].plot(data['nbytes'], data['time']/data['iterations'], marker='x', label=f'{i+1} Node{"s" if i else ""}',color = "blue",linestyle="")
    xs = np.linspace(data['nbytes'].min(), data['nbytes'].max(), 10000)
    fit = stats.linregress(data['nbytes'], data['time']/data['iterations'])
    alpha = fit.intercept
    beta = fit.slope
    axs[i].plot(xs, alpha + beta*xs,marker = "", color='red')
    axs[i].set_xlabel('Message Size (bytes)')
    axs[i].set_ylabel('Time per Ping-Pong (seconds)')
    axs[i].set_title(f'({"Same" if i == 0 else "Different"} Node), latency = %.2e s, inverse bandwidth = %.2e s/B' % (alpha, beta))
    axs[i].legend()
fig.savefig('ping_pong.png')
plt.show()
