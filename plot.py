#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 04:03:13 2018

@author: katerina

PLOTS
"""

import matplotlib.pyplot as plt
import numpy as np

'''
x, y = np.loadtxt('QW_NORM_Phiouter_qplate-Coin_n1_T1_steps25_only_NORM.txt', delimiter=',', unpack=True)
plt.plot(x, y,'ro')
plt.xlabel('iterations')
plt.ylabel('Schmidt norm')
plt.title('Schmidt norm after n=25 steps ')
plt.legend()
plt.savefig('QW_NORM_Phi_qplate-Coin_n5_T1_steps25_only_NORM.png')

'''
a, b = np.loadtxt('QW_NORM_Phiouter_qplate-Coin_n5_T1_steps25_bound_only_NORM.txt', delimiter=',', unpack=True)
plt.plot(a, b,'ro')
plt.xlabel('iterations')
plt.ylabel('Schmidt norm')
plt.title('Schmidt norm after n=25 steps  ')
plt.legend()
plt.savefig('QW_NORM_Phi_qplate-Coin_n5_T1_steps25_bound_only_NORM.png') 

