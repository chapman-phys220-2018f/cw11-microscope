#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Enea Dodi
# Student ID: 2296306
# Email: dodi@chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: cw11
###

def u(t):
    return np.array([np.cos(t),-1*np.sin(t)])

def derivSandC(x,v):
    return np.array([v,-x])

def derivu(t):
    J = np.zeros((2,2))
    J[0][0] = 0
    J[0][1] = 1
    J[1][0] = -1
    J[1][1] = 0
    return J@u(t)

def derivuTU(tk,uk):
    J = np.zeros((2,2))
    J[0][0] = 0
    J[0][1] = 1
    J[1][0] = -1
    J[1][1] = 0
    u = np.array([np.cos(tk),-1*np.sin(uk)])
    return J@u

def x(t):
    return np.cos(t)

def y(u):
    return -1*np.sin(u)

def eulers(N):
    t = np.linspace(0,10*np.pi,N)
    dt = t[1]-t[0]
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    x[0] = 1
    y[0] = 0
    for i in range(1,len(t)):
        theDerivs = derivu(t[i-1])
        x[i] = x[i-1] + dt*theDerivs[0]
        y[i] = y[i-1] + dt*theDerivs[1]
    #u=np.array([x,y])
    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, x, color = "blue")
    a.plot(t,y,color = "red")
    plt.show()
    return x

def Heun(N):
    t = np.linspace(0,10*np.pi,N)
    dt = t[1]-t[0]
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    x[0] = 1
    y[0] = 0
    for i in range(1,len(t)):
        theDerivs_two = derivu(t[i-1])
        theDerivs_one = derivu(t[i])
        x[i] = x[i-1] + (dt/2)*(theDerivs_two[0] + theDerivs_one[0])
        y[i] = y[i-1] + (dt/2)*(theDerivs_two[1] + theDerivs_one[1])
    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, x, color = "blue")
    a.plot(t,y,color = "red")
    plt.show()
    return x
        
def secondOrderRungeKutta(N):
    t = np.linspace(0,10*np.pi,N)
    dt = t[1]-t[0]
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    x[0] = 1
    y[0] = 0
    for i in range(1,len(t)):
        K1 = dt*derivu(t[i-1])
        K21 = dt*derivu(t[i-1]+dt/2)
        K22 = dt*derivu(t[i-1]+K1[1]/2)
        x[i] = (x[i-1] + K21[0])
        y[i] = (y[i-1] + K22[1])
    s = plt.figure(figsize=(8,6))
    a = plt.axes()
    a.plot(t, x, color = "blue")
    a.plot(t,y,color = "red")
    plt.show()
    return 

def fourthOrderRungeKutta(N):
    t = np.linspace(0,10*np.pi,N)
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
    a.plot(t, x, color = "blue")
    a.plot(t,y,color = "red")
    plt.show()
    return 