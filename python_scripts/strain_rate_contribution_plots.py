import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_time_strain_rate_contribution(fname):
    f = plt.figure()
    data =  pd.read_csv(fname)
    print(data)
    g=sns.lineplot(data=data, x="time",  y="dg_glide_xx", label="glide strain rate")
    g=sns.lineplot(data=data, x="time",  y="dg_climb_xx", label="climb strain rate")
    g.set(xscale='log')
    g.set(yscale='log')
    g.set(ylabel='strain rate [1/s]')
    g.set(xlabel='time [s]')
    g.legend()

if __name__ == "__main__":
    fname = "strain_rate_contribution_ph1.csv"
    plot_time_strain_rate_contribution(fname)
    plt.show()