import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
data_1_node = pd.read_csv('one-node-pp.dat', delimiter=",", header=None, names=['iterations', 'nbytes','time'])
print(data_1_node)
fig, axs = plt.subplots(1, 2, figsize=(15, 6))
axs[0].plot(data_1_node['nbytes'], data_1_node['time']/data_1_node['iterations'], marker='x', label='1 Node',color = "blue",linestyle="")
xs = np.linspace(data_1_node['nbytes'].min(), data_1_node['nbytes'].max(), 10000)
fit = stats.linregress(data_1_node['nbytes'], data_1_node['time']/data_1_node['iterations'])
alpha = fit.intercept
beta = fit.slope
print(alpha)
axs[0].plot(xs, alpha + beta*xs,marker = "", color='red')
axs[0].set_xlabel('Message Size (bytes)')
axs[0].set_ylabel('Time per Ping-Pong (seconds)')
axs[0].set_title('(Same Node), latency = %.2e s, inverse bandwidth = %.2e s/B' % (alpha, beta))
axs[0].legend()
fig.savefig('ping_pong_same_node.png')
plt.show()
