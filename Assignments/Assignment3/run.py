import os,sys
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import pandas as pd
import seaborn as sns

sns.set()

#make
out = os.popen("make").read()

dict = {
    "n_nodes": [],
    "ppn": [],  # 0 --> optimized, 1 --> standard
    "time": [],
}


data_input = pd.DataFrame.from_dict(dict)

filename = str(sys.argv[1])
outfile = "output.txt"

print(filename)

for execution in range(10):
  print("execution", execution)
  arg = "python script.py 1 2 8"
  alloc_log = os.popen(arg).read()
  for n_nodes in [1, 2]:
    for ppn in [1, 2, 4, 8]:
      arg = "mpirun -ppn " + str(ppn) + " -n " + str(ppn*n_nodes) + " -f hostfile ./a.out " + filename
      print(arg)
      out = os.popen(arg).read()
      with open(outfile, 'w') as f:
          f.write(out)
      f.close()
      print(out)
      # print(out)
      out = out.split("\n")
      out = out[2]
      data_input = data_input.append({"n_nodes": n_nodes, "ppn": ppn,  "time": out}, ignore_index=True)
  print("\n")


data_input["(n_nodes, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, data_input["n_nodes"]), map(str, data_input["ppn"])))

print(data_input)

data_input.to_csv("results.csv")

  # plot = sns.catplot(x="(P, ppn)", y="time", data=data_input[i], kind="box", col="D", hue="mode")
  # #plt.title(collectives[i], fontsize=18) #changes title of D=2048
  # plot.fig.subplots_adjust(top=0.9)
  # plot.fig.suptitle(collectives[i])
  # plt.savefig("plot_"+collectives[i]+".jpg")
  #plt.show()
os.popen("python plot.py").read()

