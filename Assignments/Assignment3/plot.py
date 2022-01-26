import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import pandas as pd
import seaborn as sns

df = pd.read_csv("results.csv")
plot = sns.catplot( x="(n_nodes, ppn)", y="time", data=df, kind="bar",  sharey=False, height=9, aspect=11/9)
plot.fig.subplots_adjust(top=0.9, bottom=0.07, left=0.07)
plot.ax.set_xlabel("(n_nodes, ppn)", fontsize=15)
plot.ax.set_ylabel("time", fontsize=15)
plt.savefig("plot.jpg")
# plt.show()
 