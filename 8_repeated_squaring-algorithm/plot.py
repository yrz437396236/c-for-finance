# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 14:42:15 2018

@author: 43739
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_excel('C:\\Users\\43739\\OneDrive\\us\\c++ coding\\8_repeated_squaring-algorithm\\5X5.xlsx',sep=',')
x1=df[['Repeated Squares_k']].values
x2=df[['Brute Force Multiplication_k']].values
x3=df[['Repeated Squares_n']].values
x4=df[['Brute Force Multiplication_n']].values
y1=df[['n']].values
y2=df[['k']].values
plt.figure()
plt.plot(y1,x1,label='Repeated Squares')
plt.plot(y1,x2,label='Brute Force Multiplication')
plt.legend()
plt.xlabel("Exponent",size=10)
plt.ylabel("Time(s)",size=10)
plt.title("Computation time vs. Exponent Plot(5X5)",size=10)
plt.savefig('C:\\Users\\43739\\OneDrive\\us\\c++ coding\\8_repeated_squaring-algorithm\\5X5_result_k.png',dpi = 3000)

plt.figure()
plt.plot(y2,x3,label='Repeated Squares')
plt.plot(y2,x4,label='Brute Force Multiplication')
plt.legend()
plt.xlabel("Size",size=10)
plt.ylabel("Time(s)",size=10)
plt.title("Computation time vs. Size of the Matrix(Exponent=100)",size=10)
plt.savefig('C:\\Users\\43739\\OneDrive\\us\\c++ coding\\8_repeated_squaring-algorithm\\5X5_result_n.png',dpi = 3000)