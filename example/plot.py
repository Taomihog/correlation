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
def linear_interpolate(x, y, x_i):
    y_i = np.interp(x_i, x, y, left = -1.0, right = -1.0)
    return y_i
x_i = [i for i in range(0,5000)]


fig,ax = plt.subplots(figsize=(16, 4))

folder_path = os.getcwd() + "/example/H5alpha"
prefix = "out"

files_with_prefix = [f for f in os.listdir(folder_path)]

# ax.plot(pyResult['# zMean'],pyResult['FMean'], "k", linewidth = 2, label = "python code")
for file in files_with_prefix:
    print(folder_path + '/' + file) 
    cppResult = pd.read_csv(folder_path + '/' + file)
    # print(cppResult.columns)
    ax.plot(cppResult['total_extension_nm'],cppResult['force_avg_pN'], 
            label = "GPU result, gene = "+file[len(prefix):-len(".csv")])


plt.ylim(9.8,18.2)
# plt.ylim(0,30)

plt.xlim(650,9000)
plt.xlabel("extension (nm)")
plt.ylabel("force (pN")
plt.legend(loc='best')
plt.savefig("result.svg", dpi = 600)
plt.show()

