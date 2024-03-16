# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 20:49:16 2023

@author: taomihog
"""
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_palette("husl")

import pandas as pd
import numpy as np
import os

plot = True
fmax = 20
fmin = 10
min_bp_trace_to_save = 1500

folder_to_import = 'H5alpha'

folder_path = os.getcwd() + "\\" + folder_to_import
files = [f for f in os.listdir(folder_path)]
# files =files[0:3]
filename_out = folder_to_import + '.txt'
open(filename_out, 'w')

def linear_interpolate(x, y, x0):
    y0 = np.interp(x0, x, y, left = -1.0, right = -1.0)
    return y0

if plot:
    fig,ax = plt.subplots(figsize=(16, 4))



for file in files:
    print(file) 
    res = pd.read_csv(folder_path + '\\' + file)
    x = res['total_extension_nm'].values
    y = res['force_avg_pN'].values
    
    # remove the high force region
    il = 0
    ih = len(x) - 1;
    imid = (int)(il + (ih - il)>>1)
    cnt = 0
    while(True):
        cnt += 1
        imid = (int)(il + (ih - il)//2)
        if y[imid] < fmax:
            il = imid
        else:
            ih = imid
        if ih - il == 1:
            break;
    
    x = x[:ih]
    y = y[:ih]
    while(ih - 2 >= 0 and y[ih - 1] - y[ih - 2] > 0):
        ih -= 1
    if ih + 20 < len(x):
        x = x[:ih + 20]
        y = y[:ih + 20]
    
    # remove the low force region
    il = 0
    ih = len(x) - 1;
    imid = (int)(il + (ih - il)>>1)
    cnt = 0
    while(True):
        cnt += 1
        imid = (int)(il + (ih - il)//2)
        if y[imid] < fmin:
            il = imid
        else:
            ih = imid
        if ih - il == 1:
            break;
    
    x = x[il:]
    y = y[il:]
    il = 0
    while(il + 1 < len(y) and y[il] - y[il + 1] < 0):
        il += 1
    if il >= 20:
        x = x[il - 20:]
        y = y[il - 20:]
    
    # interpolate the integer extension start at zero
    x -= x[0]
    xmax = x[-1]
    x0 = [x for x in range(xmax)]
    y0 = linear_interpolate(x, y, x0)
    if xmax < min_bp_trace_to_save:
        # do not save if the sequence is too short, for test purposes
        continue
    
    # check data integrity
    with open(filename_out, 'a') as file_out:
        # Write lines one by one
        
        file_out.write(file[:-len(".csv")])
        file_out.write(',')
        [file_out.write(f"{y:.4f},") for y in y0]
        file_out.write("\n")
    
    if plot:
        ax.plot(x, y, label = file[:-len(".csv")] + "_ori")
        pass

with open(filename_out, 'r') as file_in:
    Lines = file_in.readlines()
     
for cnt,line in enumerate(Lines):
    array = [float(i) for i in line.split(',')[1:-1]]
    y1 = np.array(array)
    x1 = [x for x in range(len(y1))]
    
    if plot and False:
        ax.plot(x1, y1, label = chr(cnt + 48),linestyle= ':', color = 'k')

print("count of lines: ", len(Lines))

if plot:
    # plt.ylim(5,25)
    # plt.xlim(0,9000)
    # plt.legend(loc='best')
    plt.xlabel("extension (nm)")
    plt.ylabel("force (pN")
    plt.savefig("result.svg", dpi = 600)
    plt.show()

