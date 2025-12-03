import matplotlib.pyplot as plt    
import numpy as np
import pandas as pd
import scipy.optimize as sp

def strong_model(p,f):
    """ Amdahl's law model for strong scaling
    p : number of processes
    f : fraction of serial code
    """
    return 1.0/ (f + (1-f)/p)

def weak_model(p,f):
    """ Model for weak scaling
    p : number of processes
    f : fraction of serial code
    """
    return f + (1-f)*p

def Strong():
    strong_data = pd.read_csv("pll_strong_ham.dat",delimiter=",")
    strong_data["speedup"] = strong_data["time"].iloc[0]/strong_data["time"]
    weights = 1 / strong_data["speedup"] 
    param,covar = sp.curve_fit(strong_model,strong_data["size"],strong_data["speedup"],p0=[0.9],bounds=(0,1),sigma=weights)
    fitted_f = param[0]
    xs = np.arange(1,16,0.1)
    ys = strong_model(xs,fitted_f)


    fig, ax = plt.subplots(1,1,figsize=(8,6))
    ax.plot(strong_data["size"],strong_data["speedup"],marker="x",linestyle = "")
    ax.plot(xs, ys, label=f"Fitted model (f={fitted_f:.4f})",color="red")
    #ax.set_xscale("log",base=2)
    #ax.set_yscale("log",base=2)
    ax.set_xlabel("Number of processes")
    ax.set_ylabel("Speedup")
    ax.set_title("Strong Scaling of MPI Pi Calculation (N = 10000000)")
    ax.legend()
    plt.show()
    fig.savefig("strong_scaling_fit.png")
def Weak():
    weak_data = pd.read_csv("pll_weak_ham.dat",delimiter=",")
    serial = weak_data[weak_data["size"]==1].reset_index(drop=True)
    pll = weak_data[weak_data["size"] * 1000000 == weak_data["N"]].reset_index(drop=True)
    weak_speedup = serial["time"]/pll["time"]

    xs = np.arange(1,16,0.1)
    param,covar = sp.curve_fit(weak_model,pll["size"],weak_speedup,p0=[0.2],bounds=(0,1))
    fitted_f = param[0]
    ys = weak_model(xs,fitted_f)

    fig,ax = plt.subplots(1,1,figsize=(8,6))
    ax.plot(pll["size"],weak_speedup,marker="x",linestyle = "")
    ax.plot(xs, ys, label=f"Fitted model (f={fitted_f:.4f})",color="red")
    ax.set_xlabel("Number of processes")
    ax.set_ylabel("Speedup")
    ax.set_title("Weak Scaling of MPI Pi Calculation")
    ax.legend()
    plt.show()
    fig.savefig("weak_scaling.png")
    
    efficiency = weak_speedup / pll["size"]
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    ax.plot(pll["size"],efficiency,marker="x",linestyle = "")
    ax.plot(xs, weak_model(xs,fitted_f)/xs, label=f"Fitted model (f={fitted_f:.4f})",color="red")
    ax.set_xlabel("Number of processes")
    ax.set_ylabel("Efficiency")
    ax.set_title("Weak Scaling Efficiency of MPI Pi Calculation")
    ax.legend()
    plt.show()
    fig.savefig("weak_scaling_efficiency.png")
    
    
Weak()