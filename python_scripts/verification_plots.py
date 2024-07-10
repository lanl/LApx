import numpy as np
import pandas as pd
import seaborn as sns
import os
from glob import glob
from matplotlib import pyplot as plt

def getSubDirs(currdir):
    dir_list = []
    dirs = os.listdir(".")
    for dir in dirs:
        if os.path.isdir(dir):
            if not dir.startswith("."):
                dir_list.append(dir)

    return dir_list


def getSimOptions(cdir):
    s = cdir.split("_")
    sym_type = s[0]
    stress = s[1]
    temperature = s[2]
    climb_opt = s[3]
    return sym_type, stress, temperature, climb_opt

def read_stress_strain(cdir):
    if cdir.startswith("OLD"):
        cnames = ["TIME", "EVM", "E11", "E22" ,"E33","EEVM","EE11","EE22","EE33","EPVM","EP11","EP22","EP33","DVM","D11","D22","D33","DEVM","DE11","DE22","DE33","DPVM","DP11","DP22","DP33","SVM","S11","S22","S33"]
    if cdir.startswith("NEW"):
        cnames = ["TIME", "EVM", "E11", "E22" ,"E33","EEVM","EE11","EE22","EE33","EPVM","EP11","EP22","EP33","DVM","D11","D22","D33","DEVM","DE11","DE22","DE33","DPVM","DP11","DP22","DP33","SVM","S11","S22","S33",  "avg_echVM", "avg_hyd_ech", "avg_ech11", "avg_ech22", "avg_ech33", "avg_conc", "avg_gb_conc"]

    data = pd.read_csv(cdir+'/test-str_str.out', sep = "\s+|\t+|\s+\t+|\t+\s+", header=None, skiprows=1)
    data.columns = cnames
    sym_type, stress, temperature, climb_opt = getSimOptions(cdir)
    data["symtype"] = sym_type
    data["stress"] = stress
    data["temperature"] = temperature
    data["climb_opt"] = climb_opt
    print(climb_opt)
    return data



if __name__ == "__main__":
    cwd =  os.getcwd()
    subdirs = getSubDirs(cwd)

    for i, cdir in enumerate(subdirs):
        if i == 0:
            data = read_stress_strain(cdir)
        else :
            datatemp = read_stress_strain(cdir)
            data = pd.concat([data, datatemp], axis=0, ignore_index=True)

    print(data.shape)
    idx = data["climb_opt"] == "CLIMB"
    data_no_climb = data[idx]
    print(data_no_climb.shape)

    g=sns.relplot(data=data_no_climb, kind='line',x="TIME",y="DPVM", col="temperature",row="stress", hue="symtype", style="symtype")
    g.set(yscale="log")
    g.set(xscale="log")
    g.set(xlim=[0.05,1e5])
    plt.show()
