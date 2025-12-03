import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

serial_data = pd.read_csv("serial.dat",delimiter=",")
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(serial_data['N'], serial_data["error"], label="Serial", marker="x",linestyle = "")
ax.hlines(y=0.0,xmin=serial_data['N'].min(), xmax=serial_data['N'].max(), color='red', linestyle='--', label="y =0")
ax.set_xlabel("N")
ax.set_ylabel("Error")
ax.set_title("Serial Computation Error vs N")
ax.legend()
plt.show()
fig.savefig("serial_error.png")