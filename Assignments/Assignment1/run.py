import os
import numpy as np
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Make
out = os.popen("make").read()
tsteps = 2
P_count = [16, 36, 49, 64]
N_count = [16**2, 32**2, 64**2, 128**2, 256**2, 512**2, 1024**2] 
data = np.zeros((len(P_count), len(N_count), 3, 5), np.float64) 

for exec_count in range(0, 5):
    arg = "~/UGP/allocator/src/allocator.out 64 8"
    alloc_log = os.popen(arg).read()
    #print(alloc_log)
    print("====== Execution %d ======\n" % exec_count)
    print("Allocated hosts:")
    with open("hostsimproved", 'r') as f:
        for line in f.read().split("\n"):
            print(line)
    for i,P in enumerate(P_count):
        print("Number of processes: %d" % P)
        for j, N in enumerate(N_count):
            print("N: %d" % N)
            arg = "mpirun -np " + str(P) + "  -f hostsimproved ./halo " + str(N) + " " + str(tsteps)
            output = os.popen(arg).read()
            print(output)
            output = output.split("\n")
            output = np.array([float(t) for t in output])
            data[i, j, :, exec_count] = output
        print("")
    print("")

print(data)

for i,P in enumerate(P_count):
    np.save('data' + str(P) + '.npy', data[i])

os.popen("python3 plot.py").read()
