import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import pandas as pd
import seaborn as sns

sns.set()

#make
out = os.popen("make").read()

collectives = ["Bcast", "Reduce", "Gather", "Alltoallv"]

dict = {
    "D": [],
    "P": [],
    "ppn": [],
    "mode": [],  # 0 --> optimized, 1 --> standard
    "time": [],
}

node_map = {
    "4": 1,
    "16": 2
}
data_input = []
for i in range(4):
  data_input.append(pd.DataFrame.from_dict(dict))

for execution in range(10):
  print("execution", execution)
  for P in [4, 16]:
    for ppn in [1, 8]:
      g = node_map[str(P)]
      os.popen("python script.py " + str(g) + " " + str(P/g) + " " + str(ppn)).read()
      for D in [16, 256, 2048]:
        arg = "mpirun -np " + str(P*ppn) + " -f hostfile ./coll " + str(D*128)
        out = os.popen(arg).read()
        print(out)
        #print(out)
        out = out.split("\n")
        out = np.array([float(out[i]) for i in range(len(out)-1)])
        for i in range(4):
          data_input[i] = data_input[i].append({"D": D, "P": P, "ppn": ppn, "mode": 0, "time": out[2*i]}, ignore_index=True)
          data_input[i] = data_input[i].append({"D": D, "P": P, "ppn": ppn, "mode": 1, "time": out[2*i+1]}, ignore_index=True)
  print("\n")

print(data_input)

for i in range(4):
  data_input[i]["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, data_input[i]["P"]), map(str, data_input[i]["ppn"])))

  print(data_input[i])
  data_input[i].to_csv('data_' + collectives[i] + ".csv")

  # plot = sns.catplot(x="(P, ppn)", y="time", data=data_input[i], kind="box", col="D", hue="mode")
  # #plt.title(collectives[i], fontsize=18) #changes title of D=2048
  # plot.fig.subplots_adjust(top=0.9)
  # plot.fig.suptitle(collectives[i])
  # plt.savefig("plot_"+collectives[i]+".jpg")
  #plt.show()
os.popen("python plot.py").read()
