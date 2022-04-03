import math
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

# plt.style.use('ggplot')
global coords
#coords = [[1, 1], [2, 3], [4, -1], [6, 5], [7, 0]]
# coords = [[0,4], [1,3], [2,4], [3,13],[4,36], [5, 79], [6,148],[7,249]]
coords = []
X_min = 0
X_max = 40
Y_min = 0
Y_max = 40


class BuilderPoint:
    """
    Class for GUI management
    """

    def __init__(self, graph, fig, ax, btn_vander, btn_lagrange, btn_newton_matrix, btn_newton_ddn, btn_spline, spline):
        """
        Init the class
        """
        self.release = None
        self.graph = graph
        self.fig = fig
        self.ax = ax
        self.btn_vander = btn_vander
        self.btn_lagrange = btn_lagrange
        self.btn_newton_matrix = btn_newton_matrix
        self.btn_newton_ddn = btn_newton_ddn
        self.btn_spline = btn_spline
        self.cid = None
        self.motion = None
        self.continuePlotting = False
        self.dragging_point = None
        self.spline = spline
        self.line = None

    def connect(self):
        """
        Connect the button and mouse event
        """
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.continuePlotting = True
        self.btn_vander.pack_forget()
        self.btn_lagrange.pack_forget()
        self.btn_newton_matrix.pack_forget()
        self.btn_newton_ddn.pack_forget()
        self.btn_spline.pack_forget()

    def disconnect(self):
        """
        Disconnect the button and mouse event
        :return:
        """
        self.fig.canvas.mpl_disconnect(self.cid)
        self.fig.canvas.mpl_disconnect(self.motion)
        self.fig.canvas.mpl_disconnect(self.release)
        self.continuePlotting = False
        if coords:
            self.btn_vander.pack()
            self.btn_lagrange.pack()
            self.btn_newton_matrix.pack()
            self.btn_newton_ddn.pack()
            self.btn_spline.pack()

    def drag_and_drop(self):
        """
        Connect the mouse event
        """
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)

    def change_state(self):
        """
        If spline is active, enable picking
        If plotting is active, disable it
        else, enable Plotting
        """
        if self.spline.upd:
            self.drag_and_drop()
        elif self.continuePlotting:
            self.disconnect()
        else:
            self.connect()

    def onclick(self, event):
        """ callback method for mouse click event
        :type event: MouseEvent
        """
        x, y = event.xdata, event.ydata
        if self.spline.upd:
            self.dragging_point = self.find_neighbor_point(event)
        else:
            x = round(x)
            y = round(y)
            coords.append([x, y])
        self.update_plot()

    def on_motion(self, event):
        """ callback method for mouse motion event
        :type event: MouseEvent
        """
        if not self.dragging_point:
            return
        if event.xdata is None or event.ydata is None:
            return
        coords[:] = [[round(event.xdata), round(event.ydata)] if e == self.dragging_point else e for e in coords]
        self.dragging_point = [round(event.xdata), round(event.ydata)]
        self.update_plot()

    def on_release(self, event):
        """ callback method for mouse release event
        :type event: MouseEvent
        """
        if self.dragging_point:
            self.dragging_point = None
            self.update_plot()

    def find_neighbor_point(self, event):
        """ Find point around mouse position
        :rtype: ((int, int)|None)
        :return: (x, y) if there are any point around mouse else None
        """
        distance_threshold = 1.5
        nearest_point = None
        min_distance = math.sqrt(2 * (100 ** 2))
        for x, y in coords:
            distance = math.hypot(event.xdata - x, event.ydata - y)
            if distance < min_distance:
                min_distance = distance
                nearest_point = [x, y]
        if min_distance < distance_threshold:
            return nearest_point
        return None

    def update_plot(self):
        """
        Update the current plot
        """
        # Update current plot
        if not self.line:
            self.line = self.ax.scatter(coords[0][0], coords[0][1], c='b', s=50)
        else:
            self.line.set_offsets(coords)
        self.graph.draw()
        if self.spline.upd:
            self.spline.show_matrix(True)


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
        Construct the Newton matrix
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
                print("i", i)
                print(self.newton_matrix_div[j] - self.newton_matrix_div[j - 1])
                print(self.matrix_of_x[j] - self.matrix_of_x[j - i - 1])
                self.newton_matrix_div[j] = (self.newton_matrix_div[j] - self.newton_matrix_div[j - 1]) / (
                        self.matrix_of_x[j] - self.matrix_of_x[j - i - 1])
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


class Spline:
    """
    Class for Spline interpolation
    """

    def __init__(self, root, graph, ax):
        self.root = root
        self.graph = graph
        self.ax = ax
        self.matrix_of_x = None
        self.matrix_of_y = None
        self.spline = None
        self.coef = None
        self.show = False
        self.upd = False
        self.curve = None

    def matrix_x_and_y(self):
        """
        Divide the coords into two separate list of coords x and y
        :rtype: [],[]
        :return: list of x and list of y
        """

        self.matrix_of_x, self.matrix_of_y = build_matrix_x_and_y(False)

    def build_spline(self):
        """
        Construct the BSpline curve
        :return: An array of values representing the spline function
        """
        # Interpolate Spline parameter : Find the B-spline representation of an N-D curve.
        # Modify k to change the degree
        tck, u = interpolate.splprep([self.matrix_of_x, self.matrix_of_y], k=3, s=0)

        # Get more precision
        u = np.linspace(0, 1, num=100, endpoint=True)

        # Evaluate a B-spline
        self.coef = interpolate.splev(u, tck)

    def show_matrix(self, drag=False):
        """
        Manage the GUI of the Bspline interpolation
        :return: the interpolation of a BSpline
        """

        if not self.show or drag:
            self.show = True
            self.upd = True
            self.matrix_x_and_y()
            self.build_spline()

            # Print Poly
            if not self.curve:
                self.curve, = self.ax.plot(self.coef[0], self.coef[1], 'cyan')
            else:
                self.curve.set_data(self.coef[0], self.coef[1])
            self.graph.draw()

        else:
            self.show = False
            self.upd = False
            self.curve = None


def horner(x_matrix, coef, x):
    # Initialize result
    result = 0

    # Evaluate value of polynomial
    # using Horner's method
    for i in range(len(x_matrix) - 1, -1, -1):
        result = result * (x - x_matrix[i]) + coef[i]
    return result


def build_matrix_x_and_y(spline=True):
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

    if not spline:
        list_x = np.append(list_x, list_x[0])
        list_y = np.append(list_y, list_y[0])
        return list_x, list_y
    return np.asarray(list_x), np.asarray(list_y)


def sort_coords():
    """
    Sort the coords list with the x coordonnates
    :return: Sorted list
    """
    coords.sort(key=lambda coord: coord[0])


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
    spline = Spline(root, graph, ax)

    btn_vander = Button(root, text="Watch Vander Matrix", command=vandermonde.show_matrix, bg="blue", fg="white")
    btn_lagrange = Button(root, text="Print graph lagrange", command=lagrange.show_matrix, bg="green", fg="white")
    btn_newton_matrix = Button(root, text="Print graph newton matrix", command=newton.show_matrix, bg="purple",
                               fg="white")
    btn_newton_ddn = Button(root, text="Print graph newton ddn", command=newton.show_matrix_ddn, bg="magenta",
                            fg="white")
    btn_spline = Button(root, text="Print spline cubique", command=spline.show_matrix, bg="cyan",
                        fg="white")
    builder_point = BuilderPoint(graph, fig, ax, btn_vander, btn_lagrange, btn_newton_matrix, btn_newton_ddn,
                                 btn_spline, spline)

    # Manage button for builder point
    b = Button(root, text="Start/Stop or picking", command=builder_point.change_state, bg="red", fg="white")
    b.pack()

    root.mainloop()


if __name__ == '__main__':
    init()
