import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
from scipy.optimize import curve_fit
import os
script_dir = os.path.dirname(os.path.abspath(__file__))




#Note, we cannot just plot (p,t(p)) to estimate f given that t(1) changes as we increase p given that we also increase n with p. 

def straight(x,a,b):
    return a + b*x

#import data
weak_data = pd.read_csv("/home/ciaran/HPC/hello-world/weak.dat",delimiter=",",names = ["p","t_p","supplied_t_s"] )
P = np.array(weak_data["p"])
T = np.array(weak_data["t_p"])
T_s = np.array(weak_data["supplied_t_s"])
epsilon = np.linspace(0,10e-2,num = len(T))
gen_T_s = 0.4*T[0] + epsilon*P
betas = T_s/T
f_average = np.mean(betas)
#print(np.mean(betas))
#print(gen_T_s)
#generate fits
params_1, covar_1 = curve_fit(straight,T,T_s,[0,1])
f = params_1[1]
intercept = params_1[0]
f_err = np.sqrt(np.diag(covar_1))[1]

#generate analyticals
xs = np.arange(0,3,0.1)
ys_fit = straight(xs,intercept,f)
xs_p = np.arange(1,256,0.5)
ys_average = straight(xs,0,f_average)
#speedup
speedup_data =  (T_s + P*(T-T_s))/T
speedup_pred = f_average+(1-f_average)*xs_p
#plot data
fig, axs = plt.subplots(1,2,figsize = (10,10))
axs[0].plot(T,T_s,color = "blue",label = "weak data",marker ="X",linestyle = "")
axs[0].plot(xs,ys_fit,color = "red",label = "fit (using curve fit): beta(p,p) = "+str(round(f,2)))
axs[0].plot(xs,ys_average,color ="green",label = "fit using average: beta(p,p) = "+str(round(f_average,2)))
axs[0].set_xlabel("Total Runtime, T(p,p)")
axs[0].set_ylabel("Serial Runtime, Ts(p)/P")
axs[0].legend()
#plot axis 2
axs[1].plot(P,speedup_data,color = "blue",label = "speedup - data",linestyle = "",marker = "X")
axs[1].plot(xs_p,speedup_pred,color = "red",label = "speedup - fit")
axs[1].set_xlabel("N Processors, p")
axs[1].set_ylabel("Speedup")
axs[1].legend()
save_path = os.path.join(script_dir, "weak_fit.png")
fig.savefig(save_path)
plt.show()