from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

# plt.style.use('ggplot')
global coords
#coords = [[1, 1], [2, 3], [4, -1], [6, 5], [7, 0]]
#coords = [[0,4], [1,3], [2,4], [3,13],[4,36], [5, 79], [6,148],[7,249]]
coords = []
X_min = 0
X_max = 40
Y_min = 0
Y_max = 40


class BuilderPoint:
    """
    Class for GUI management
    """

    def __init__(self, graph, fig, ax, btn_vander, btn_lagrange, btn_newton_matrix, btn_newton_ddn):
        self.graph = graph
        self.fig = fig
        self.ax = ax
        self.btn_vander = btn_vander
        self.btn_lagrange = btn_lagrange
        self.btn_newton_matrix = btn_newton_matrix
        self.btn_newton_ddn = btn_newton_ddn
        self.cid = None
        self.continuePlotting = False

    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.continuePlotting = True
        self.btn_vander.pack_forget()
        self.btn_lagrange.pack_forget()
        self.btn_newton_matrix.pack_forget()
        self.btn_newton_ddn.pack_forget()

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid)
        self.continuePlotting = False
        if coords:
            self.btn_vander.pack()
            self.btn_lagrange.pack()
            self.btn_newton_matrix.pack()
            self.btn_newton_ddn.pack()

    def change_state(self):
        if self.continuePlotting:
            self.disconnect()
        else:
            self.connect()

    def onclick(self, event):
        x, y = event.xdata, event.ydata
        x = round(x)
        y = round(y)
        coords.append([x, y])

        self.ax.scatter([x], [y], c='b', s=50)
        self.graph.draw()


class Vandermonde:
    """
    Class for Vandermonde interpolation
    """

    def __init__(self, root, graph, ax):
        self.root = root
        self.graph = graph
        self.ax = ax
        self.matrix_of_x = None
        self.matrix_of_y = None
        self.vander = None
        self.coef = None
        self.inverse_vander = None
        self.show = False
        self.label = None

    def matrix_x_and_y(self):
        """
        Divide the coords into two separate list of coords x and y
        :return: list of x and list of y
        """

        self.matrix_of_x, self.matrix_of_y = build_matrix_x_and_y()

    def build_vandermonde(self):
        """
        Construct the Vandermonde matrix
        :return: matrix
        """
        self.vander = np.vander(self.matrix_of_x, increasing=True)

    def calcul_coefficient(self):
        """
        Compute the coef of the interpolation
        :return: matrix of coef
        """
        try:
            # Inverse the vandermonde matrix
            self.inverse_vander = np.linalg.inv(self.vander)
            # Compute the coef
            self.coef = self.inverse_vander.dot(self.matrix_of_y)
        except ValueError as err:
            print("Nombre trop grand", err)

    def show_matrix(self):
        """
        Manage the GUI of the vandermonde interpolation and calculate the coef with the vandermonde methods
        :return: the interpolation
        """
        if not self.show:
            self.show = True
            sort_coords()
            self.matrix_x_and_y()
            self.build_vandermonde()
            self.calcul_coefficient()

            # Print result in the label
            self.label = LabelFrame(self.root, text="Vandermonde matrix")
            self.label.pack(fill="both", expand="yes")
            Label(self.label, text=self.vander).pack(fill=BOTH)
            Label(self.label, text="Coef :\n").pack()
            Label(self.label, text=self.coef).pack(fill=BOTH)

            # Get poly
            x = np.linspace(X_min, X_max, 100)
            y = np.array([np.sum(np.array([self.coef[i] * (j ** i) for i in range(len(self.coef))])) for j in x])

            # Print Poly
            self.ax.plot(x, y, linewidth=2.0, c="blue")
            self.graph.draw()

        else:
            self.show = False
            self.label.pack_forget()


class Lagrange:
    """
    Class for lagrange interpolation
    """

    def __init__(self, root, graph, ax):
        self.root = root
        self.graph = graph
        self.ax = ax
        self.matrix_of_x = None
        self.matrix_of_y = None
        self.coef = 0
        self.show = False

    def matrix_x_and_y(self):
        """
        Divide the coords into two separate list of coords x and y
        :return: list of x and list of y
        """

        self.matrix_of_x, self.matrix_of_y = build_matrix_x_and_y()

    def Lagrange_n_k(self):
        """
        Calculate the coefs with lagrange methods
        :return: Coefs of lagrange
        """
        # Create a one polynomial class
        X = np.poly1d([1, 0])
        self.coef = 0

        # Compute the coefs
        for n in range(len(self.matrix_of_x)):
            product_coef = 1
            for k in range(len(self.matrix_of_x)):
                if k != n:
                    product_coef *= (X - self.matrix_of_x[k]) / (self.matrix_of_x[n] - self.matrix_of_x[k])
            self.coef += product_coef * self.matrix_of_y[n]

    def show_matrix(self):
        """
        Manage the GUI of the lagrange interpolation and calculate the coef with the lagrange methods
        :return: the interpolation of lagrange
        """
        if not self.show:
            sort_coords()
            self.show = True
            self.matrix_x_and_y()
            self.Lagrange_n_k()

            # Print Poly
            x = np.linspace(X_min, X_max, 100)
            self.ax.plot(x, np.polyval(self.coef, x), linewidth=2.0, c="green")
            self.graph.draw()
        else:
            self.show = False


class Newton:

    def __init__(self, root, graph, ax):
        self.root = root
        self.graph = graph
        self.ax = ax
        self.matrix_of_x = None
        self.matrix_of_y = None
        self.coef = 0
        self.show = False
        self.show_div = False
        self.newton_matrix = None
        self.newton_matrix_div = []

    def matrix_x_and_y(self):
        """
        Divide the coords into two separate list of coords x and y
        :return: list of x and list of y
        """

        self.matrix_of_x, self.matrix_of_y = build_matrix_x_and_y()

    def build_newton_matrix(self):
        """
        Construct the Vandermonde matrix
        :return: matrix
        """
        self.newton_matrix = np.zeros((len(self.matrix_of_x), len(self.matrix_of_x)), dtype=int)
        for i in range(len(self.matrix_of_x)):
            for j in range(len(self.matrix_of_x)):
                if j == 0:
                    self.newton_matrix[i][j] = 1
                elif j > i:
                    break
                else:
                    self.newton_matrix[i][j] = self.newton_matrix[i][j - 1] * (
                                self.matrix_of_x[i] - self.matrix_of_x[j - 1])

    def calcul_coefficient_matrix(self):
        """
        Compute the coef of the interpolation
        :return: matrix of coef
        """
        try:
            # Inverse the Newton matrix
            inverse_newton = np.linalg.inv(self.newton_matrix)
            # Compute the coef
            self.coef = inverse_newton.dot(self.matrix_of_y)

        except ValueError as err:
            print("Nombre trop grand", err)

    def diff_div(self):
        print(self.matrix_of_y)
        for i in range(len(self.matrix_of_y)):
            self.newton_matrix_div.append(self.matrix_of_y[i])

        for i in range(len(self.matrix_of_x)):
            for j in range(len(self.matrix_of_x) - 1, i, -1):
                print("i" , i)
                print(self.newton_matrix_div[j] - self.newton_matrix_div[j - 1])
                print(self.matrix_of_x[j] - self.matrix_of_x[j - i - 1])
                self.newton_matrix_div[j] = (self.newton_matrix_div[j] - self.newton_matrix_div[j - 1]) / (self.matrix_of_x[j] - self.matrix_of_x[j - i - 1])
                print("after result")
                print(self.newton_matrix_div[j])
        print(self.newton_matrix_div)


    def show_matrix(self):
        if not self.show:
            self.show = True
            sort_coords()
            self.matrix_x_and_y()
            self.build_newton_matrix()
            self.calcul_coefficient_matrix()
            print(self.coef)

            # Get poly
            x = np.linspace(X_min, X_max, 100)
            y = []
            for i in x:
                y.append(horner(self.matrix_of_x, self.coef, i))

            # Print Poly
            self.ax.plot(x, y, linewidth=2.0, c="purple")
            self.graph.draw()

        else:
            self.show = False

    def show_matrix_ddn(self):
        if not self.show_div:
            self.show_div = True
            sort_coords()
            self.matrix_x_and_y()
            self.diff_div()

            # Get poly
            x = np.linspace(X_min, X_max, 100)
            y = []
            for i in x:
                y.append(horner(self.matrix_of_x, self.newton_matrix_div, i))

            # Print Poly
            self.ax.plot(x, y, linewidth=2.0, c="magenta")
            self.graph.draw()

        else:
            self.show_div = False


def horner(x_matrix, coef, x):
    # Initialize result
    result = 0

    # Evaluate value of polynomial
    # using Horner's method
    for i in range(len(x_matrix) - 1, -1, -1):
        result = result * (x - x_matrix[i]) + coef[i]
    return result


def build_matrix_x_and_y():
    """
    Divide the coords into two separate list of coords x and y
    :return: list of x and list of y
    """

    # Create list of X coordinates from coords
    list_x = []
    for i in range(len(coords)):
        list_x.append(coords[i][0])

    # Create a list of Y coordinates from coords
    list_y = []
    for i in range(len(coords)):
        list_y.append(coords[i][1])
    return np.asarray(list_x), np.asarray(list_y)


def sort_coords():
    """
    Sort the coords list with the x coordonnates
    :return: Sorted list
    """
    sorted(coords, key=lambda coord: coords[0])


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
    ax.set_xlim(X_min, X_max)
    ax.set_ylim(Y_min, Y_max)
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()

    # Add the figure of mathplotlib to the tkinter interface
    graph = FigureCanvasTkAgg(fig, master=root)
    graph.get_tk_widget().pack(side="top", fill='both', expand=True)

    # Init Builder Point
    vandermonde = Vandermonde(root, graph, ax)
    lagrange = Lagrange(root, graph, ax)
    newton = Newton(root, graph, ax)

    btn_vander = Button(root, text="Watch Vander Matrix", command=vandermonde.show_matrix, bg="blue", fg="white")
    btn_lagrange = Button(root, text="Print graph lagrange", command=lagrange.show_matrix, bg="green", fg="white")
    btn_newton_matrix = Button(root, text="Print graph newton matrix", command=newton.show_matrix, bg="purple",
                               fg="white")
    btn_newton_ddn = Button(root, text="Print graph newton ddn", command=newton.show_matrix_ddn, bg="magenta",
                            fg="white")
    builder_point = BuilderPoint(graph, fig, ax, btn_vander, btn_lagrange, btn_newton_matrix, btn_newton_ddn)

    # Manage button for builder point

    b = Button(root, text="Start/Stop", command=builder_point.change_state, bg="red", fg="white")
    b.pack()

    root.mainloop()


if __name__ == '__main__':
    init()
