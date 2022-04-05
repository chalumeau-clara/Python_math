import math
from tkinter import Tk, Label

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import numpy.linalg
from scipy.misc import derivative
from scipy.stats import expon
import sympy as sp


def function(x):
    return x * x * x + x - 1 # expon.ppf(-x)


def print_function(a, b):
    list_x = numpy.linspace(a, b, 100)
    plt.plot(list_x, function(list_x), "k-")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


def variation(f, a, b):
    '''
    True => Croissant
    False => Decroissant
    '''
    list_x = np.linspace(a, b, 100)
    d = derivative(f, list_x, dx=1e-6)
    if d[0] < 0:
        for i in d:
            if i > 0:
                print("Non strict sur [" + str(a) + ',' + str(b) + ']')
                return
        print("Décroissant sur [" + str(a) + ',' + str(b) + ']')
        return False
    if d[0] > 0:
        for i in d:
            if i < 0:
                print("Non strict sur [" + str(a) + ',' + str(b) + ']')
                return
        print("Croissant sur [" + str(a) + ',' + str(b) + ']')
        return True


def convexe_or_concave(f, a, b):
    '''
    True => convexe
    False => concave
    '''
    list_x = np.linspace(a, b, 100)
    d = derivative(f, list_x, dx=1e-6, n=2)
    if d[0] < 0:
        for i in d:
            if i > 0:
                print("Non strict sur [" + str(a) + ',' + str(b) + ']')
                return
        print("Concave sur [" + str(a) + ',' + str(b) + ']')
        return False
    if d[0] > 0:
        for i in d:
            if i < 0:
                print("Non strict sur [" + str(a) + ',' + str(b) + ']')
                return
        print("Convexe sur [" + str(a) + ',' + str(b) + ']')
        return True


def lagrange_method(f, a, b):
    list_x = numpy.linspace(a, b, 100)
    plt.plot(list_x, function(list_x), "k-")

    if f(a) * derivative(f, a, n=2) > 0:
        pivot = a
        x0 = b
    else:
        pivot = b
        x0 = a
    count, xn = 10, x0
    while count > 0:
        x0 = xn
        xn = x0 - f(x0) * (x0 - pivot) / (f(x0) - f(pivot))
        plt.plot(xn, 0, linestyle='none', marker='o', c='lime',
                    markersize=5)
        count -= 1
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
    return xn


def newton_method(f, a, b):
    list_x = numpy.linspace(a, b, 100)
    plt.plot(list_x, function(list_x), "k-")

    if f(a) * derivative(f, a, n=2) > 0:
        x0 = a
    else:
        x0 = b
    count, xn = 10, x0
    while count > 0:
        x0 = xn
        xn = x0 - f(x0) / derivative(f, x0)
        plt.plot(xn, 0, linestyle='none', marker='o', c='lime',
                 markersize=5)
        count -= 1
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
    return xn


if __name__ == '__main__':
    a, b = 0, 2
    print_function(a, b)
    var = variation(function, a, b)
    if var != True and var != False:
        print("La fonction n'est pas monotone sur cet interval")
        exit()
    conv = convexe_or_concave(function, a, b)
    if conv != True and conv != False:
        print("La fonctione n'est pas strictement convexe ou concave sur cet interval")
        exit()
    print("Methode par Lagrange, x = " + str(lagrange_method(function, a, b)))
    print("Methode par Newton, x = " + str(newton_method(function, a, b)))
