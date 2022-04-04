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
xcoords = [0, 1, 2, 3, 4, 6, 7, 9, 11, 12, 15, 16, 17, 18, 20]
ycoords = [2, 1, 0, -1, -3, -1, 0, 2, 4, 5, 7, 10, 8, -3, -10]

class Approximation:

    def __init__(self, ax, graph):
        self.nb_point = 15
        self.epsilon = 0
        self.degree = 1
        self.sum_xi = []
        self.sum_xc = []
        self.matrix = []
        self.min_x = 0
        self.max_x = 0
        self.ax = ax
        self.graph = graph

    def init_point(self, ax):
        # User input
        self.nb_point = 15 #int(input("Please, enter the number of points:\n"))
        self.epsilon = 2 #float(input("Enter the incertitude:\n"))
        self.min_x = 0 #int(input("Enter the min axis for x:\n"))
        self.max_x = 22 #int(input("Enter the max axis for x:\n"))
        min_y = -12 #int(input("Enter the min axis for y:\n"))
        max_y = 12 #int(input("Enter the max axis for y:\n"))

        # Set axis
        ax.set_xlim(self.min_x - 3, self.max_x + 3)
        ax.set_ylim(min_y - 3, max_y + 3)

        # Build Lists
        for i in range(self.nb_point):
            xcoords.append(random.uniform(self.min_x, self.max_x))
        xcoords.sort()

        ycoords.append(random.uniform(min_y, max_y))
        for i in range(self.nb_point - 1):
            rand = random.uniform(ycoords[-1] - self.epsilon, ycoords[-1] + self.epsilon)
            if rand < min_y: rand = min_y
            if rand > max_y: rand = max_y
            ycoords.append(rand)
        print("Coordonnates of x ", xcoords)
        print("Coordonnates of y ", ycoords)

    def init_ex(self, ax):
        # User input
        self.nb_point = 15  # int(input("Please, enter the number of points:\n"))
        self.epsilon = 2  # float(input("Enter the incertitude:\n"))
        self.min_x = 0  # int(input("Enter the min axis for x:\n"))
        self.max_x = 52  # int(input("Enter the max axis for x:\n"))
        min_y = -62  # int(input("Enter the min axis for y:\n"))
        max_y = 12  # int(input("Enter the max axis for y:\n"))

        # Set axis
        ax.set_xlim(self.min_x - 3, self.max_x + 3)
        ax.set_ylim(min_y - 3, max_y + 3)

    def get_degree(self):
        print("Le degree d'un polynome est égal au nombre de variation")
        print("Coordonnates of x ", xcoords)
        print("Coordonnates of y ", ycoords)
        for i in range(self.nb_point - 2):
            if ycoords[i] < ycoords[i + 1] and ycoords[i + 1] > ycoords[i + 2] \
                    or ycoords[i] > ycoords[i + 1] and ycoords[i + 1] < ycoords[i + 2]:
                self.degree += 1
        print(self.degree)

    def sum_power_of_x(self):
        for i in range(2 * self.degree, -1, -1):
            self.sum_xi.append(self.somme(i))

    def somme(self, degree):
        sum = 0
        for i in xcoords:
            sum += math.pow(i, degree)
        return sum

    def sum_of_power_xc(self):
        for i in range(self.degree, -1, -1):
            self.sum_xc.append(self.somme_xc(i))

    def somme_xc(self, degree):
        sum = 0
        for i in range(len(xcoords)):
            sum += math.pow(xcoords[i], degree) * ycoords[i]
        return sum

    def build_matrix(self):
        for i in range(self.degree + 1):
            ligne = []
            for j in range(i, i + self.degree + 1):
                ligne.append(self.sum_xi[j])
            self.matrix.append(ligne)

    def show_matrix(self):
        self.get_degree()
        self.sum_power_of_x()
        self.sum_of_power_xc()
        self.build_matrix()
        print(self.matrix)
        print(self.sum_xc)
        A = np.array(self.matrix)
        B = np.array(self.sum_xc)
        X = np.linalg.solve(A, B)
        print(X)

        x = 22
        y = 25
        z = 50

        rx = np.polyval(X, x)
        ry = np.polyval(X, y)
        rz = np.polyval(X, z)
        print("P(22) = ", rx, " P(25) = ", ry, "P(50) = ", rz)

        x = np.linspace(self.min_x, self.max_x, 100)
        self.ax.plot(x, np.polyval(X, x), linewidth=2.0, c="green")
        self.graph.draw()

    def show_matrix_ex(self):
        self.get_degree()
        self.sum_power_of_x()
        self.sum_of_power_xc()
        self.build_matrix()
        print(self.matrix)
        print(self.sum_xc)
        A = np.array(self.matrix)
        B = np.array(self.sum_xc)
        X = np.linalg.solve(A, B)
        print(X)

        x = 22
        y = 25
        z = 50

        rx = np.polyval(X, x)
        ry = np.polyval(X, y)
        rz = np.polyval(X, z)
        print("P(22) = ", rx, " P(25) = ", ry, "P(50) = ", rz)


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

    # Add the figure of mathplotlib to the tkinter interface
    graph = FigureCanvasTkAgg(fig, master=root)
    graph.get_tk_widget().pack(side="top", fill='both', expand=True)

    # Add a plot
    ax = fig.add_subplot(111)
    app = Approximation(ax, graph)
    app.init_ex(ax)
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()


    ax.plot(xcoords, ycoords, '.')
    print("Coordonnates of x ", xcoords)
    print("Coordonnates of y ", ycoords)
    graph.draw()

    app.show_matrix()

    root.mainloop()


if __name__ == '__main__':
    init()
