# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 14:42:15 2018

@author: 43739
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_excel('C:\\Users\\43739\\OneDrive\\us\\c++ coding\\7_simulation exercises\\0.18_0.2_result.xlsx',sep=',')
x1=df[['theoretical']].values
x2=df[['simulated']].values
y=df[['n']].values
plt.figure()
plt.plot(y,x1,label='Theory')
plt.plot(y,x2,label='Experiment')
plt.legend()
plt.xlabel("Number of Coin Tosses in Each Game",size=10)
plt.ylabel("Alice’s Probability of Winning",size=10)
plt.title("q=0.18 p=0.2",size=10)
plt.savefig('C:\\Users\\43739\\OneDrive\\us\\c++ coding\\7_simulation exercises\\0.18_0.2_result.png',dpi = 3000)