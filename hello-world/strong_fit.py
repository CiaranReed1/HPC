import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
from scipy.optimize import curve_fit
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def amdahls(p,t1,f): #define scaling model
    return t1*f + t1*(1-f)/p

#import data
strong_data = pd.read_csv("/home/ciaran/HPC/hello-world/strong.dat",delimiter=",",names = ["p","t_p"] )
ps = np.array(strong_data["p"])
ts = np.array(strong_data["t_p"])

#scipy curve fit
params, covar = curve_fit(amdahls,ps,ts,[64,1])
t1 = params[0]
f = params[1]
errors = np.sqrt(np.diag(covar))
ana_xs = np.arange(1,260,0.5)
ana_ys = amdahls(ana_xs,t1,f) 

#speedup
speedup_data = ts[0] / ts
speedup_pred = 1/ (f + ((1-f)/ana_xs))
max_speedup = 1/f
#plot - axis 1 
fig, axs = plt.subplots(1,2,figsize = (10,10))
axs[0].plot(ps,ts,color = "blue",marker = "X",linestyle = "",label = "strong data")
axs[0].plot(ana_xs,ana_ys,color = "red",label = "fit: t1 = "+str(round(t1,1))+", f = "+str(round(f,3)))
axs[0].set_xlabel("Number of processors, p")
axs[0].set_ylabel("Runtime")
axs[0].set_title("Runtime against p")
axs[0].legend()

#plot - axis 2
axs[1].plot(ps,speedup_data,color = "blue",marker = "X",linestyle = "",label = "strong data")
axs[1].plot(ana_xs,speedup_pred,color = "red",label = "predicted speedup")
axs[1].set_xlabel("Number of processors, p")
axs[1].set_ylabel("Speedup")
axs[1].set_title("speedup against p, max predicted speedup = "+str(round(max_speedup,2)))
axs[1].legend()
#save and show
save_path = os.path.join(script_dir, "strong_fit.png")
fig.savefig(save_path)
plt.show()