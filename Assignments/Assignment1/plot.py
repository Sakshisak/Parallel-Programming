import os
import numpy as np
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P_count = [16, 36, 49, 64]

for P in P_count:
  send1 = []
  packed1 = []
  type1 = []
  data = np.load('data' + str(P) + '.npy')

  for k in range(7):
    d_send = [t for t in data[k, 0, :]]
    d_packed = [t for t in data[k, 1, :]]
    d_type = [t for t in data[k, 2, :]]
    send1.append(d_send)
    packed1.append(d_packed)
    type1.append(d_type)

  mean_send = np.sum(np.array(send1), axis=1)/5
  mean_packed = np.sum(np.array(packed1), axis=1)/5
  mean_type = np.sum(np.array(type1), axis=1)/5
  # print("send1: \n",send1)
  # print("packed1: \n", packed1)
  # print("type1: \n", type1)

  #for k in range(7):
   # send1[k] = [np.log10(t) for t in send1[k]]
   # packed1[k] = [np.log10(t) for t in packed1[k]]
   # type1[k] = [np.log10(t) for t in type1[k]]

  fig_data = []
  fig_data.append(send1)
  fig_data.append(packed1)
  fig_data.append(type1)

  fig = plt.figure(figsize = (10,8))
  ax = fig.add_subplot(111)
  bplot_send = ax.boxplot(send1, vert=True, showfliers=False, patch_artist=True, boxprops = dict(color='black', facecolor = 'CornflowerBlue'), positions=np.array(range(7))*4.0 + 1)
  bplot_packed = ax.boxplot(packed1, vert=True, showfliers=False, patch_artist=True, boxprops= dict(color='black', facecolor='DarkSalmon'), positions=np.array(range(7))*4.0 + 2)
  bplot_type = ax.boxplot(type1, vert=True, showfliers=False, patch_artist=True, boxprops= dict(color='black', facecolor='yellow'), positions=np.array(range(7))*4.0 + 3)
  #boxprops = dict(linestyle='--', linewidth=2, color='Black', facecolor = 'red', alpha = .4)
  # colors = ['#0000FF', '#00FF00', '#FFFF00']
  line_send = ax.plot(np.array(range(7))*4.0 + 1, mean_send, color="blue")
  line_packed = ax.plot(np.array(range(7))*4.0 + 2, mean_packed, color="red")
  line_type = ax.plot(np.array(range(7))*4.0 + 3, mean_type, color="orange")
  for median in bplot_send['medians']: 
    median.set(color ='blue', 
               linewidth = 1)
  
  for median in bplot_packed['medians']: 
    median.set(color ='red', 
               linewidth = 1)
               
  for median in bplot_type['medians']: 
    median.set(color ='orange', 
               linewidth = 1)

  # x-axis labels 
  ax.set_xticklabels(['16^2', '32^2',  
                    '64^2', '128^2', '256^2', '512^2', '1024^2'])
  #ax.set_yticklabels(['%.2f'%(10**t) for t in ax.get_yticks()])
  # axes labels
  ax.set_xticks(np.array(range(7))*4.0 + 2)
  ax.set_xlabel('No. of data entries', fontsize = 14)           
  ax.set_ylabel('Time in seconds', fontsize = 14)

  ax.legend([bplot_send["boxes"][0], bplot_packed["boxes"][0], bplot_type["boxes"][0], line_send[0], line_packed[0], line_type[0]], ["halo_send", "halo_packed", "halo_type", "halo_send", "halo_packed", "halo_type"], loc='upper left') 
  plt.title("Process count : "+ str(P), fontsize=18)
  
 # ax.get_xaxis().tick_bottom() 
 # ax.get_yaxis().tick_left()

  plt.show()
  plt.savefig("plot"+str(P)+".jpg")
