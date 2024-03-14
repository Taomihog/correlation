# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 04:05:56 2024

@author: xiang
"""

import numpy as np

def linear_interpolate(x, y, x_i):
    y_i = np.interp(x_i, x, y, left = -1.0, right = -1.0)
    return y_i

# Example usage
x = [0.0, 1.0, 2.0, 3.0, 4.0]
y = [1.0, 2.0, 3.0, 4.0, 5.0]
x_i = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]

# Interpolate y values for x_i
y_i = linear_interpolate(x, y, x_i)

# Print interpolated y values
print(y_i)