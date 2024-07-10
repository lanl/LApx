import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_stress_strain(fname):
    # Plot the stress strain Curve
    f = plt.figure()
    data =  pd.read_csv(fname)
    g=sns.lineplot(data=data, x="tot_strain_xx",  y="stress_xx")
    g.set(xlabel='strain [mm/mm]')
    g.set(ylabel='stress_xx [MPa]')

def plot_stress_time(fname):
    # Plot the stress strain Curve
    f = plt.figure()
    data =  pd.read_csv(fname)
    g=sns.lineplot(data=data, x="time",  y="stress_xx")
    g.set(xlabel='time [s]')
    g.set(ylabel='stress_xx [MPa]')

def plot_time_strain_rate(fname):
    # Plot strain rate versus time
    f = plt.figure()
    data =  pd.read_csv(fname)
    g=sns.lineplot(data=data, x="time",  y="tot_strain_rate_xx", label="total strain rate")
    g=sns.lineplot(data=data, x="time",  y="elastic_strain_rate_xx", label="elastic strain rate")
    g=sns.lineplot(data=data, x="time",  y="plastic_strain_rate_xx", label="plastic strain rate")
    g.set(xscale='log')
    g.set(yscale='log')
    g.set(ylabel='strain rate [1/s]')
    g.set(xlabel='time [s]')
    g.legend()



if __name__ == "__main__":
    fname = "simulation_macro_averages.csv"
    plot_stress_strain(fname)
    plot_time_strain_rate(fname)
    plot_stress_time(fname)
    plt.show()