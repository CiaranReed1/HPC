import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

serial_data = pd.read_csv("serial.dat",delimiter=",")
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(np.log(serial_data['N']),np.log(serial_data["error"]), label="Serial", marker="x",linestyle = "")
res = stats.linregress(np.log(serial_data['N']), np.log(serial_data['error']))
xs = np.linspace(np.log(serial_data['N']).min(), np.log(serial_data['N']).max(), 1000)
ax.plot(xs, res.intercept + res.slope*xs, color='red', label="Fit: log(Error) = %.2f + %.2f * log(N)" % (res.intercept, res.slope))
#ax.hlines(y=0.0,xmin=np.log(serial_data['N']).min(), xmax=np.log(serial_data['N']).max(), color='red', linestyle='--', label="y =0")
ax.set_xlabel("log(N)")
ax.set_ylabel("log(Error)")
ax.set_title("Serial Computation Error vs N")
ax.legend()
plt.show()
fig.savefig("serial_error.png")

fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(serial_data['N'], serial_data["time"], label="Serial", marker="x",linestyle = "")
res = stats.linregress(serial_data['N'], serial_data['time'])
xs = np.linspace(serial_data['N'].min(), serial_data['N'].max(), 1000)
ax.plot(xs, res.intercept + res.slope*xs, color='red', label="Fit: time = %.2e + %.2e * N" % (res.intercept, res.slope))
ax.set_xlabel("N")
ax.set_ylabel("Time (s)")
ax.set_title("Serial Computation Time vs N")
ax.legend()
plt.show()
fig.savefig("serial_time.png")
