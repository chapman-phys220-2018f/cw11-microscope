#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Enea DOdi
# Student ID: 2296306
# Email: dodi@chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: c209
###

def gradient(x):
        """gradient function
        Args: x - input numpy array consiting of a set of n domain points.
        Returns: a n by n matrix creating a finite difference linear operator"""
        dx = x[1]-x[0] #dx should be a very small number
        ones = np.ones(x.size-1) #an array contiang just ones, which will be the values placed in the gradient
        diag1 = np.diag(-1*ones,-1) #The negative one which will represent (fi - 1)
        diag2 = np.diag(ones,1)#Postive integer one which will represent (fi - 1)
        gradient = diag1+diag2
        
        #The parts of the gradient where central difference cannot be used
        gradient[0][0] = -2 
        gradient[0][1] = 2
        gradient[-1][-2] = -2
        gradient[-1][-1] = 2
        
        return np.divide(gradient,2*dx)

def deriv(x,f):
    """deriv function
    Args: x - input numpy array consisting of a set of n domain points
          f - function which the derivative will be calculated from
          Returns: DerivY - the y values of the derivative"""
    func = np.vectorize(f)
    yVals = func(x) #Get y values of the domain points after function f is applied to them
    #derivative creation!
    D = gradient(x)
    DerivY = D@yVals
    #Return the y values of the derivative.
    return DerivY

def deriv_vs_normal_plot(x,f,name):
        """deriv_vs_normal_plot(x,f)
        Args: x - input numpy array consisting of a set of n domain points
        f - function which the derivtive will be calculated from
        Returns: nothing but returns 3 graphs, one which is simply the function f(x), 
        the other is simply the derivative f'(x), and the third is a graph with both the 
        plotting of the functino f(x) and f'(x)"""
        Df = deriv(x,f)
        titl = name+"(x)"
        titl_deriv = name + "'(x)"
        func = np.vectorize(f)
        yVal = func(x)
        s = plt.figure(figsize=(8,6))
        a = plt.axes()
        a.plot(x,yVal, color = "blue")
        a.set(xlabel= "X values", ylabel= "Y values", title=titl)
        plt.show()
        s = plt.figure(figsize=(8,6))
        a = plt.axes()
        a.plot(x,Df, color = "red")
        a.set(xlabel= "X values", ylabel= "Y values", title=titl_deriv)
        plt.show()
        s = plt.figure(figsize=(8,6))
        a = plt.axes()
        a.plot(x,yVal, color = "blue")
        a.plot(x,Df,color = "red")
        a.set(xlabel= "X values", ylabel= "Y values", title=titl + " vs " +titl_deriv)
        plt.show()
        