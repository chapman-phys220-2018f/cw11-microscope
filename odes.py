#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Enea Dodi
# Student ID: 2296306
# Email: dodi@chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: Cw11
###
'''Contains all the code for cw11'''

import numpy as np
from matplotlib import pyplot as plt
import array_calculus as ac

def u(t):
    '''u function
    Args - t time value
    returns derivative at point t
    '''
    return np.array([np.cos(t),-1*np.sin(t)])

def derivu(t):
    '''derivu function
    Args - t time value
    returns the derivative of the values at t'''
    J = np.zeros((2,2))
    J[0][0] = 0
    J[0][1] = 1
    J[1][0] = -1
    J[1][1] = 0
    return J@u(t)

def x(t):
    '''x function
    args - t time value
    returns cos at t
    '''
    return np.cos(t)

def y(u):
    '''y function
    args - t time value
    returns -sin at t
    '''
    return -1*np.sin(u)

def deriv(x,v):
    '''deriv function
    args - x the cosine
           v the -sin
    returns the derivatives of the cos and -sin valuess
    '''
    J = np.zeros((2,2))
    J[0][0] = 0
    J[0][1] = 1
    J[1][0] = -1
    J[1][1] = 0
    u = np.array([x,v])
    return J@u

def Eulers(N,r0):
    '''eulers function
    Args: N - the number of time steps
          r0 - the inital values 
    returns the graphs using the points calculated by euler formula
    for x(t) and v(t)
    '''
    t = np.arange(0,10*(np.pi),2*(np.pi)/N)
    r = np.zeros([2,len(t)])
    dt = t[1]-t[0]

    r[0][0] = r0[0]
    r[1][0] = r0[1]
    for i in range(1,len(t)):
        theDerivs = derivu(t[i-1])
        r[0][i] = r[0][i-1] + dt*theDerivs[0]
        r[1][i] = r[1][i-1] + dt*theDerivs[1]
    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, r[0], color = "blue", label = "Cosine")
    a.plot(t,r[1],color = "red", label = "Sine")
    a.legend(loc='upper left')
    plt.show()


def Heun(N,r0):
    '''Heun function
    Args: N - the number of time steps
          r0 - the inital values 
    returns the graphs using the points calculated by the Heun formula
    for x(t) and v(t)
    '''
    t = np.arange(0,10*(np.pi),2*(np.pi)/N)
    r = np.zeros([2,len(t)])
    dt = t[1]-t[0]

    r[0][0] = r0[0]
    r[1][0] = r0[1]
    for i in range(1,len(t)):
        HeunOne = derivu(t[i-1])
        proxX = r[0][i-1] + (dt)*(HeunOne[0])
        proxV = r[1][i-1] +(dt)*(HeunOne[1])
        u = deriv(proxX,proxV)
        r[0][i] = r[0][i-1] + (dt/2)*(HeunOne[0]+u[0])
        r[1][i] = r[1][i-1] + (dt/2)*(HeunOne[1]+u[1])

    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, r[0], color = "blue", label = "Cosine")
    a.plot(t,r[1],color = "red", label = "Sine")
    a.legend(loc='upper left')
    plt.show()
    
def SecondOrderRungeKutta(N,r0):
    '''SecondOrderRungeKutta function
    Args: N - the number of time steps
          r0 - the inital values 
    returns the graphs using the points calculated by the SecondOrderRungeKutta formula
    for x(t) and v(t)
    '''
    t = np.arange(0,10*(np.pi),2*(np.pi)/N)
    dt = t[1]-t[0]
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    x[0] = r0[0]
    y[0] = r0[1]
    for i in range(1,len(t)):
        K1 = dt*derivu(t[i-1])
        K21 = dt*derivu(t[i-1]+dt/2)
        K22 = dt*derivu(t[i-1]+K1[1]/2)
        x[i] = (x[i-1] + K21[0])
        y[i] = (y[i-1] + K22[1])
    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, x, color = "blue", label = "Cosine")
    a.plot(t,y,color = "red", label = "Sine")
    a.legend(loc='upper left')
    plt.show()

def FourthOrderRungeKutta(N,r0):
    '''FourthOrderRungeKutta function
    Args: N - the number of time steps
          r0 - the inital values 
    returns the graphs using the points calculated by the FourthOrderRungeKutta formula
    for x(t) and v(t)
    '''
    t = np.arange(0,10*(np.pi),2*(np.pi)/N)
    dt = t[1]-t[0]
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    x[0] = 1
    y[0] = 0
    for i in range(1,len(t)):
        K1 = dt*derivu(t[i-1])
        K21 = dt*derivu(t[i-1]+dt/2)
        K22 = dt*derivu(t[i-1]+K1[1]/2)
        K32 = dt*derivu(t[i-1]+K22[1]/2)
        K41 = dt*derivu(t[i-1]+dt)
        K42 = dt*derivu(t[i-1]+K32[1])
        x[i] = x[i-1] + (K1[0]+2*K21[0]+2*K21[0]+K41[0])/6
        y[i] = y[i-1] + (K1[1]+2*K22[1]+2*K32[1]+K42[1])/6
    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, x, color = "blue", label = "Cosine")
    a.plot(t,y,color = "red", label = "Sine")
    a.legend(loc='upper left')
    plt.show()