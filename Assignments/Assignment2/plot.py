import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import pandas as pd
import seaborn as sns

collectives = ["Bcast", "Reduce", "Gather", "Alltoallv"]
for i in range(4):
    df = pd.read_csv("data_" + collectives[i] + ".csv")
    plot = sns.catplot( x="(P, ppn)", y="time", data=df, kind="bar", col="D", hue="mode", sharey=False)
    plot.fig.subplots_adjust(top=0.9)
    plot.fig.suptitle(collectives[i])
    plt.savefig("plot_"+collectives[i]+".jpg")
    # plt.show()