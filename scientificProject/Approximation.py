import math
import random
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

global xcoords
global ycoords

xcoords = []
ycoords = []


def init_point(ax):
    # User input
    nb_point = int(input("Please, enter the number of points:\n"))
    epsilon = float(input("Enter the incertitude:\n"))
    min_x = int(input("Enter the min axis for x:\n"))
    max_x = int(input("Enter the max axis for x:\n"))
    min_y = int(input("Enter the min axis for y:\n"))
    max_y = int(input("Enter the max axis for y:\n"))

    # Set axis
    ax.set_xlim(min_x - 3, max_x + 3)
    ax.set_ylim(min_y - 3, max_y + 3)

    # Build Lists
    for i in range(nb_point):
        xcoords.append(random.uniform(min_x, max_x))
    xcoords.sort()

    ycoords.append(random.uniform(min_y, max_y))
    for i in range(nb_point - 1):
        rand = random.uniform(ycoords[-1] - epsilon, ycoords[-1] + epsilon)
        if rand < min_y: rand = min_y
        if rand > max_y: rand = max_y
        ycoords.append(rand)
    print("Coordonnates of x ", xcoords)
    print("Coordonnates of y ", ycoords)


def init():
    """
    Initialize the GUI
    :return:
    """
    # Initialise a window, its background and its dimension.
    root = Tk()
    root.config(background='white')
    root.geometry("1000x700")

    # Give title to the window
    Label(root, text="Interpolation", bg='white').pack()

    # Init the matplotlib space
    fig = Figure()

    # Add a plot
    ax = fig.add_subplot(111)
    init_point(ax)
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()

    # Add the figure of mathplotlib to the tkinter interface
    graph = FigureCanvasTkAgg(fig, master=root)
    graph.get_tk_widget().pack(side="top", fill='both', expand=True)

    ax.plot(xcoords, ycoords, '.')
    graph.draw()
    # Init Builder Point

    # Manage button for builder point

    root.mainloop()


if __name__ == '__main__':
    init()
