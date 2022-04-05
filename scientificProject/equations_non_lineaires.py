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
    return x * x * x + x - 1
    # x * x * x + x - 1
    # 1 - x * x
    # x*x*x - 2*x + 3
    # x - np.exp(-x)
    # np.power(x,1/3)


def print_function(a, b):
    list_x = numpy.linspace(a, b, 100)
    plt.plot(list_x, function(list_x), "k-")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
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
    count, count_lines, xn = 10, 0, x0
    while count > 0:
        x0 = xn
        xn = x0 - f(x0) * (x0 - pivot) / (f(x0) - f(pivot))
        if count_lines < 3:
            plt.plot([pivot,x0], [f(pivot), f(x0)], linestyle='-', c='lime',
                    linewidth=1)
            count_lines += 1
        count -= 1
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.show()
    return xn


def newton_method(f, a, b):
    list_x = numpy.linspace(a, b, 100)
    plt.plot(list_x, function(list_x), "k-")

    if f(a) * derivative(f, a, n=2) > 0:
        x0 = a
    else:
        x0 = b
    count, count_lines, xn = 10, 0, x0
    while count > 0:
        x0 = xn
        xn = x0 - f(x0) / derivative(f, x0)
        if count_lines < 3:
            plt.plot([x0, xn], [f(x0), 0], linestyle='-', c='lime',
                 linewidth=1)
            count_lines += 1
        count -= 1
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.show()
    return xn

def is_passing_0(f, a, b):
    if f(a) * f(b) < 0:
        return True
    return False

def decoupage(f, a, b):
    final_a = a
    final_b = b
    while not is_passing_0(f, final_a, final_b) and (final_a < b and final_b > a):
        final_b -= 0.5
    if final_b <= final_a:
        final_a = a
        final_b = b
    while not is_passing_0(f, final_a, final_b) and (final_a < b and final_b > a):
        final_a += 0.5
    if final_a >= final_b:
        final_a = a
        final_b = b
    while not is_passing_0(f, final_a, final_b) and (final_a < b and final_b > a):
        final_a += 0.5
        final_b -= 0.5
    if final_a >= b:
        return False
    while is_passing_0(f, final_a, final_b) and (final_a < b and final_b > a):
        var = variation(f, final_a, final_b)
        if var == True or var == False:
            conv = convexe_or_concave(f, final_a, final_b)
            if conv == True or conv == False:
                return (final_a, final_b)
        final_b -= 0.5
    final_a = a
    final_b = b
    while is_passing_0(f, final_a, final_b) and (final_a < b and final_b > a):
        var = variation(f, final_a, final_b)
        if var == True or var == False:
            conv = convexe_or_concave(f, final_a, final_b)
            if conv == True or conv == False:
                return (final_a, final_b)
        final_a += 0.5
    final_a = a
    final_b = b
    while is_passing_0(f, final_a, final_b) and (final_a < b and final_b > a):
        var = variation(f, final_a, final_b)
        if var == True or var == False:
            conv = convexe_or_concave(f, final_a, final_b)
            if conv == True or conv == False:
                return (final_a, final_b)
        final_a += 0.5
        final_b += 0.5

if __name__ == '__main__':
    a, b = -2, 2
    print_function(a, b)
    var = variation(function, a, b)
    if var != True and var != False:
        print("La fonction n'est pas monotone sur l'intervalle ["+str(a)
              +","+str(b)+']')
        if decoupage(function, a, b) != False:
            a, b = decoupage(function,a,b)
        else:
            print("Cette fonction ne passe par 0 sur ["+str(a)
              +","+str(b)+']')
            exit()
    conv = convexe_or_concave(function, a, b)
    if conv != True and conv != False:
        print("La fonction n'est pas strictement convexe ou "
              "concave sur l'interalle ["+str(a)
              +","+str(b)+']')
        a, b = decoupage(function, a, b)
    print("Methode par Lagrange, x = " + str(lagrange_method(function, a, b)))
    print("Methode par Newton, x = " + str(newton_method(function, a, b)))
