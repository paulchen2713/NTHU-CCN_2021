# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 16:14:25 2022
@author: Paul
@file:   selection_combining.py
@dependencies: 
    conda env:  pytest
    python:     3.7.11
    numpy:      1.21.2
    matplotlib: 3.5.0
"""
import numpy as np
import matplotlib.pyplot as plt
import os

gamma_ratio_dB = np.arange(start=-10, stop=41, step=2)
Ns = [1, 2, 3] # number of received signal paths

gamma_ratio = 10**(gamma_ratio_dB/10) #Average SNR/SNR threshold in dB

plt.figure(1)
for N in Ns:
    P_outage = (1 - np.exp(-1 / gamma_ratio))**N
    plt.semilogy(gamma_ratio_dB, P_outage, label='N=' + str(N))

plt.xlim(-10, 40)
plt.ylim(0.000001, 1.5)
plt.title('theorical Selection combining outage probability')
plt.xlabel(r'$10log_{10}\left(\Gamma/\gamma_t\right)$')
plt.ylabel('Outage probability')
plt.grid(True)

# first check whether the specified path is an existing file
path = './Selection_Combining/SC.png'
# if no file exist, then save current file
if os.path.isfile(path) is False:
    plt.savefig('./Selection_Combining/SC.png')

plt.show()
