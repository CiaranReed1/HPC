import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('weak_timing.dat', sep='\t', header=None, names=['threads', 'time', 'N', 'adotb'])
data["speedup"] = (data["time"].iloc[0] / data["time"]) 
f_est = 1/1000001
speedup_pred = f_est+(1-f_est)*data["threads"]
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(data['threads'], data['time'], marker='o',color = "blue")
#ax.plot(data['threads'], speedup_pred, marker='x', color='red', linestyle='--', label='Predicted Speedup')
ax.set_xlabel('Number of Threads')
ax.set_ylabel('Speedup')
ax.set_title('Weak Scaling Performance')
fig.savefig('weak_scaling.png')
plt.show()