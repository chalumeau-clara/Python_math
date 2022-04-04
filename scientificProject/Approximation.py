import math
import random
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

global xcoords
global ycoords

xcoords = []
ycoords = []


class Approximation:
    """
    Approximation class
    :return Moindre carrée approximation
    """

    def __init__(self, ax, graph, root):
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
        self.x = 0
        self.root = root

    def init_point(self, ax):
        """
        Init the paremeters randomly
        :param ax: ax for plot
        """
        # User input
        self.nb_point = int(input("Please, enter the number of points:\n"))
        self.epsilon = float(input("Enter the incertitude:\n"))
        self.min_x = int(input("Enter the min axis for x:\n"))
        self.max_x = int(input("Enter the max axis for x:\n"))
        min_y = int(input("Enter the min axis for y:\n"))
        max_y = int(input("Enter the max axis for y:\n"))

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

    def init_ex(self, ax):
        """
        Init the paremeters by N in the exercice
        :param ax: ax for plot
        """
        global xcoords
        global ycoords
        xcoords = [0, 1, 2, 3, 4, 6, 7, 9, 11, 12, 15, 16, 17, 18, 20]
        ycoords = [2, 1, 0, -1, -3, -1, 0, 2, 4, 5, 7, 10, 8, -3, -10]

        # User input
        self.nb_point = 15
        self.epsilon = 2
        self.min_x = 0
        self.max_x = 52
        min_y = -62
        max_y = 12
        self.x = 1

        # Set axis
        ax.set_xlim(self.min_x - 3, self.max_x + 3)
        ax.set_ylim(min_y - 3, max_y + 3)

    def init_ex7(self, ax):
        """
        Init the paremeters according to ex 7
        :param ax: ax for plot
        """
        global xcoords
        global ycoords
        xcoords = [1, 2, 4, 6, 7]
        ycoords = [1, 3, 5, 10, 15]

        # User input
        self.nb_point = 5
        self.epsilon = 2
        self.min_x = 0
        self.max_x = 10
        min_y = 0
        max_y = 20

        # Set axis
        ax.set_xlim(self.min_x - 3, self.max_x + 3)
        ax.set_ylim(min_y - 3, max_y + 3)

    def init_ex8(self, ax):
        """
         Init the paremeters according to ex 8
         :param ax: ax for plot
         """
        global xcoords
        global ycoords
        xcoords = [0, 1, 2, 4, 6]
        ycoords = [7, 4, 2, 4, 12]

        # User input
        self.nb_point = 5
        self.epsilon = 2
        self.min_x = 0
        self.max_x = 10
        min_y = 0
        max_y = 15

        # Set axis
        ax.set_xlim(self.min_x - 3, self.max_x + 3)
        ax.set_ylim(min_y - 3, max_y + 3)

    def get_degree(self):
        """
        Calculate the degree of the polynome
        :return: The degree of the polynome
        """
        print("Le degree d'un polynome est égal au nombre de variation")
        for i in range(self.nb_point - 2):
            if ycoords[i] < ycoords[i + 1] and ycoords[i + 1] > ycoords[i + 2] \
                    or ycoords[i] > ycoords[i + 1] and ycoords[i + 1] < ycoords[i + 2]:
                self.degree += 1
        #print(self.degree)

    def sum_power_of_x(self):
        """
        Compute sum of xi
        :rtype: [float]
        :return: sum of x for each power
        """
        for i in range(2 * self.degree, -1, -1):
            self.sum_xi.append(self.somme(i))

    def somme(self, degree):
        """
        Compute sum of xi ** degree
        :param degree: the degree of power
        :rtype: float
        :return: sum of xi ** degree
        """
        sum = 0
        for i in xcoords:
            sum += math.pow(i, degree)
        return sum

    def sum_of_power_xc(self):
        """
        Compute sum of xc (where c is the y coord)
        :rtype: [float]
        :return: sum of xc for each power
        """
        for i in range(self.degree, -1, -1):
            self.sum_xc.append(self.somme_xc(i))

    def somme_xc(self, degree):
        """
         Compute sum of xc ** degree
         :param degree: the degree of power
         :rtype: float
         :return: sum of xc ** degree
         """
        sum = 0
        for i in range(len(xcoords)):
            sum += math.pow(xcoords[i], degree) * ycoords[i]
        return sum

    def build_matrix(self):
        """
        Build the moindre carre matrix
        :rtype: [[float]]
        :return: moindre carre matrix
        """
        for i in range(self.degree + 1):
            ligne = []
            for j in range(i, i + self.degree + 1):
                ligne.append(self.sum_xi[j])
            self.matrix.append(ligne)

    def show_matrix(self):
        """
        Show the resulting graph
        """
        self.get_degree()
        self.sum_power_of_x()
        self.sum_of_power_xc()
        self.build_matrix()
        #print(self.matrix)
        #print(self.sum_xc)
        A = np.array(self.matrix)
        B = np.array(self.sum_xc)
        X = np.linalg.solve(A, B)
        print("coef:")
        print(X)

        if self.x == 1:
            x = 22
            y = 25
            z = 50

            rx = np.polyval(X, x)
            ry = np.polyval(X, y)
            rz = np.polyval(X, z)

            # Print result in the label
            self.label = LabelFrame(self.root, text="Prévision")
            self.label.pack(fill="both", expand="yes")
            Label(self.label, text="P(22) = ").pack(fill=BOTH)
            Label(self.label, text=rx).pack(fill=BOTH)
            Label(self.label, text="P(25) = ").pack(fill=BOTH)
            Label(self.label, text=ry).pack(fill=BOTH)
            Label(self.label, text="P(50) = ").pack(fill=BOTH)
            Label(self.label, text=rz).pack(fill=BOTH)

        x = np.linspace(self.min_x, self.max_x, 100)
        self.ax.plot(x, np.polyval(X, x), linewidth=2.0, c="green")
        self.graph.draw()


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
    app = Approximation(ax, graph, root)
    user = int(input("Enter 1 for generating an approximation(by default), 2 for N define, 3 for Ex 7 and 4 for Ex 8"))
    if user == 4:
        app.init_ex8(ax)
    elif user == 2:
        app.init_ex(ax)
    elif user == 3:
        app.init_ex7(ax)
    else:
        app.init_point(ax)
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()

    # Show point
    ax.plot(xcoords, ycoords, '.')
    graph.draw()

    # Show graph
    app.show_matrix()

    root.mainloop()


if __name__ == '__main__':
    init()
