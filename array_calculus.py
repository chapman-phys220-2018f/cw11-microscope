'''
array_calculus:
    Contains functions for determining derivatives of mathematical functions using
    linear algebra.'''

import math
from matplotlib import pyplot as plt
import numpy as np

def do_inner(x, d):
    '''
    do_inner(x, d):
        Fills inner values of matrice with appropriate values. Not for
        use outside of gradient(x, f) function.
        
        Args:
            x: One dimensional numpy array of length n
            d: Two dimensional numpy array of zeros with size (n x n)
        Returns:
            d: Zeros for row k where k = [1,n-1] replaced with -1 and 1 at
            column indices n-1 and n+1 respectively
    
    '''
    n = 1
    
    while n < d.shape[0]-1:
        d[n][n-1] = -1
        d[n][n+1] = 1
        n+=1
        
    return d

def gradient(x):
    '''
    gradient(x):
        Takes 1-dimensional array x of length n, returns central
        difference matrix from x
        
        Args:
            x: One dimensional numpy array of length n containing x values
        Returns:
            d: Two dimensional numpy array of size n by n containing
            central difference (and edge case) values
    '''
    n = x.size
    d = np.zeros((n, n))
    dx = x[1]-x[0]
    d = do_inner(x,d)
    d[0][0], d[0][1] = -2, 2
    d[n-1][n-2], d[n-1][n-1] = -2, 2
    d = (d/(2*dx))
    return d

def gen_y(f, x):
    '''
    gen_y(f, x):
        Takes 1-dimensional array x and any function
        f and returns a vector with function output
        values stored in the diagonal
        
        Args:
            f: Function to evaluate
            x: Domain for function f
        Returns:
            y: n by n numpy array of f(x) values stored diagonally
    '''
    f = np.vectorize(f) # in case function isn't compatible with np arrays
    y = f(x)
    
    return y

def plot_derivatives(x, f):
    '''Helper method to display graphs of f and f'(x)
    
        Args: 
            x (np array): x values for graph
            f (function): function to evaluate
        Returns:
            None'''
    
    # f(x)
    plt.figure(figsize=(8,6))
    plt.title("f(x)")
    plt.plot(x,f(x))
    plt.show()
    
    # f'(x)
    gradient_values = gradient(x)
    y_values = gen_y(f, x)
    deriv_values = gradient_values@y_values
    
    plt.figure(figsize=(8,6))
    plt.plot(x, deriv_values)
    plt.title('f\'(x)')
    plt.show