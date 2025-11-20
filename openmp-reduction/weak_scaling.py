import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('weak_timing.dat', sep='\t', header=None, names=['threads', 'time', 'N', 'adotb']) 
f_est = 9/10
data["efficiency"] = (data['time'].iloc[0] / data['time'])  * (1/data["threads"])
data["measured_speedup"] = (data['time'].iloc[0] * data["threads"])/ data['time']
data["speedup_pred"] = f_est + (1- f_est)*data['threads']
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(data['threads'], data['measured_speedup'], marker='x',color = "blue", linestyle='', label='Data Speedup')
ax.plot(data['threads'], data['speedup_pred'], color='red', label='Gustafson-Barsis Speedup')
ax.set_xlabel('Number of Threads')
ax.set_ylabel('Speedup')
ax.legend()
ax.set_title('Weak Scaling Performance, estimated f = '+str(f_est))
fig.savefig('weak_scaling.png')
plt.show()