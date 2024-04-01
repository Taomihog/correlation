# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 23:51:38 2024

@author: xiang
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read the CSV file into a DataFrame
df = pd.read_csv('corr.csv', header=None)

# Convert the DataFrame to a NumPy array
data_array = df.values


# Create a heatmap
plt.imshow(data_array, cmap='hot', interpolation='nearest')

# Add color bar
plt.colorbar()