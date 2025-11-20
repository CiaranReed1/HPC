import pandas as pd
import matplotlib.pyplot as plt

data_files=  ['weak_timing_first_try.dat', 'weak_timing_spread.dat', 'weak_timing_close.dat']
title_suffixes = ['default', 'spread', 'close']
f_est = [9/10, 9/10,9/10]
fig, axs = plt.subplots(1,3,figsize=(8,6))
for i,file in enumerate(data_files):
    data = pd.read_csv(file, sep='\t', header=None, names=['threads', 'time', 'N', 'adotb']) 
    data["efficiency"] = (data['time'].iloc[0] / data['time'])  * (1/data["threads"])
    data["measured_speedup"] = (data['time'].iloc[0] * data["threads"])/ data['time']
    data["speedup_pred"] = f_est[data_files.index(file)] + (1- f_est[data_files.index(file)])*data['threads']
    axs[i].plot(data['threads'], data['speedup_pred'], color='red', label='Gustafson-Barsis Speedup')
    axs[i].plot(data['threads'], data['measured_speedup'], marker='x', linestyle='', color='blue')
    axs[i].set_xlabel('Number of Threads')
    axs[i].set_ylabel('Speedup')
    title_suffix = title_suffixes[i]
    axs[i].set_title('Weak Scaling Efficiency ('+title_suffix+' binding), f = '+str(f_est[i]))

    axs[i].legend()
fig.savefig('weak_scaling_speedup'+'.png')
plt.show()
