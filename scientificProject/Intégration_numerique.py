from tkinter import *
from matplotlib.widgets import TextBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

X_min = 0
X_max = 40
Y_min = 0
Y_max = 40

# Global variables for numerical integration
global integration_interval, a_text, b_text, n_text, function, numerical_integration
integration_input = [0, 0, 0]
function = "sin"

global results
results = {'rectangle_inferieur': 0,
           'rectangle_superieur': 0,
           'trapeze': 0,
           'simpson': 0}


def apply_function(func_name, x_coords):
    """
    Sets function to use for graph.
    :return: Returns application of function on x values given.
    """
    # TODO: use button(s) instead to switch between functions
    switcher = {
        "sin": np.sin,
        "cos": np.cos,
        "exp": np.exp,
    }

    if func_name == "1/sin":
        x = [1 / el for el in x_coords]
        return list(map(switcher.get("sin"), x))

    res = list(map(switcher.get(func_name), x_coords))

    return res


class NumericalIntegration:
    """
    Class for numerical integration
    """

    def __init__(self, root, graph, ax):
        self.root = root
        self.graph = graph
        self.ax = ax

    def set_a(self):
        if a_text:
            integration_input[0] = int(a_text.get(0.0, END))
            print("Set: " + str(integration_input[0]))

    def set_b(self):
        if b_text:
            integration_input[1] = int(b_text.get(0.0, END))
            print("Set: " + str(integration_input[1]))

    def set_nb_points(self):
        if n_text:
            integration_input[2] = int(n_text.get(0.0, END)) - 1
            print("Set: " + str(integration_input[2]))

    def plot_function(self):
        """
        Uses input data as well as function set by user to draw graph.
        :return:
        """
        # Clear axes
        self.ax.clear()

        x = np.linspace(integration_input[0], integration_input[1], integration_input[2])

        # Gets function from switcher if exists else uses default function
        y = np.sin(x)
        results['rectangle_inferieur'] = numerical_integration.rectangles_inferieurs("sin", integration_input[0],
                                                                                     integration_input[1],
                                                                                     integration_input[2],
                                                                                     None)
        results['rectangle_superieur'] = numerical_integration.rectangles_superieurs("sin", integration_input[0],
                                                                                     integration_input[1],
                                                                                     integration_input[2],
                                                                                     None)

        results['simpson'] = numerical_integration.simpson("sin", integration_input[0],
                                                           integration_input[1],
                                                           integration_input[2])

        results['trapeze'] = numerical_integration.rectangles_trapeze("sin", integration_input[0],
                                                                      integration_input[1],
                                                                      integration_input[2],
                                                                      None)

        self.display_results()

    def plot_last_question(self):
        """
        Uses input data as well as function set by user to draw graph.
        :return:
        """
        # Clear axes
        self.ax.clear()

        x = np.linspace(integration_input[0], integration_input[1], integration_input[2])

        # Gets function from switcher if exists else uses default function
        y = np.sin(1 / x)
        results['rectangle_inferieur'] = numerical_integration.rectangles_inferieurs("1/sin", integration_input[0],
                                                                                     integration_input[1],
                                                                                     integration_input[2],
                                                                                     None)
        results['rectangle_superieur'] = numerical_integration.rectangles_superieurs("1/sin", integration_input[0],
                                                                                     integration_input[1],
                                                                                     integration_input[2],
                                                                                     None)
        results['trapeze'] = numerical_integration.rectangles_trapeze("1/sin", integration_input[0],
                                                                      integration_input[1],
                                                                      integration_input[2],
                                                                      None)

        results['simpson'] = numerical_integration.simpson("1/sin", integration_input[0],
                                                           integration_input[1],
                                                           integration_input[2])
        # Readjust plot axis
        # self.ax.autoscale()

        # self.ax.plot(x, y, c="blue")
        # self.graph.draw()
        self.display_results()

    def rectangles_inferieurs(self, func_name, a, b, n, y_coords):
        # Clear axes
        self.ax.clear()

        h = (b - a) / n
        x_coords = list((a + (i * h) for i in range(n)))
        mapping = y_coords if y_coords else apply_function(func_name, x_coords)

        # Readjust plot axis and draw graph
        self.ax.autoscale()

        self.ax.scatter(x_coords, mapping, c="green")
        self.graph.draw()

        return h * sum(mapping)

    def rectangles_superieurs(self, func_name, a, b, n, y_coords):
        # Clear axes
        self.ax.clear()

        h = (b - a) / n
        x_coords = list((a + (i * h) for i in range(1, n + 1)))
        mapping = y_coords if y_coords else apply_function(func_name, x_coords)

        # Readjust plot axis and draw graph
        self.ax.autoscale()

        self.ax.scatter(x_coords, mapping, c="green")
        self.ax.plot(x_coords, mapping)
        self.graph.draw()

        return h * sum(mapping)

    def rectangles_trapeze(self, func_name, a, b, n, y_coords):
        # Clear axes
        self.ax.clear()

        h = (b - a) / n
        x_coords = list((a + (i * h) for i in range(1, n)))
        mapping = y_coords[1:n] if y_coords else apply_function(func_name, x_coords)

        # Readjust plot axis and draw graph
        self.ax.autoscale()

        self.ax.scatter(x_coords, mapping, c="green")
        self.ax.plot(x_coords, mapping)

        self.graph.draw()

        return ((h / 2) * (y_coords[0] + y_coords[n] + (2 * sum(mapping)))) \
            if y_coords else ((h / 2) * (apply_function(func_name, [a])[0] + apply_function(func_name, [b])[0] + (
                2 * sum(mapping))))

    def simpson_with_mapping(self, a, b, n, y_coords):
        # Clear axes
        self.ax.clear()

        if n % 2 != 0:
            print("The method can only be applied on an odd number of points")
            return
        h = (b - a) / n
        x_coords = list((a + (i * h) for i in range(n + 1)))

        mapping1 = (y_coords[i] for i in range(1, n, 2))
        mapping2 = (y_coords[i] for i in range(2, n, 2))

        # Readjust plot axis and draw graph
        self.ax.autoscale()

        self.ax.scatter(x_coords, y_coords, c="green")
        self.ax.plot(x_coords, y_coords)

        self.graph.draw()

        return (h / 3) * (y_coords[0] + y_coords[n] + (4 * sum(mapping1)) + (2 * sum(mapping2)))

    def simpson(self, func_name, a, b, n):
        # Clear axes
        self.ax.clear()

        if n % 2 != 0:
            print("The method can only be applied on an odd number of points")
            return
        h = (b - a) / n
        x_coords = list((a + (i * h) for i in range(n + 1)))

        odd_index_points = (x_coords[i] for i in range(1, n, 2))
        even_index_points = (x_coords[i] for i in range(2, n, 2))

        mapping1 = apply_function(func_name, odd_index_points)
        mapping2 = apply_function(func_name, even_index_points)

        # Readjust plot axis and draw graph
        self.ax.autoscale()

        mapping3 = apply_function(func_name, x_coords)
        self.ax.scatter(x_coords, mapping3, c="green")
        self.ax.plot(x_coords, mapping3)

        self.graph.draw()

        return (h / 3) * (apply_function(func_name, [a])[0] + apply_function(func_name, [x_coords[n]])[0] +
                          (4 * sum(mapping1)) + (2 * sum(mapping2)))

    def display_results(self):
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr = "\nRec inférieur: " + str(results['rectangle_inferieur']) + \
                  "\nRec supérieur: " + str(results['rectangle_superieur']) + \
                  "\nTrapeze: ", str(results['trapeze']) + "\nSimpson: " + str(results['simpson'])

        textstr2 = '\n'.join((
            r'Rectangle inférieur: ', str(results['rectangle_inferieur']),
            r'Rectangle supérieur: ', str(results['rectangle_superieur']),
            r'Trapeze: ', str(results['trapeze']),
            r'Simpson: ', str(results['simpson'])))

        print("\nRec inférieur: " + str(results['rectangle_inferieur']) + \
              "\nRec supérieur: " + str(results['rectangle_superieur']) + \
              "\nTrapeze: ", str(results['trapeze']) + "\nSimpson: " + str(results['simpson']))

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        self.ax.text(0.05, 0.95, textstr2, transform=self.ax.transAxes, fontsize=14,
                     verticalalignment='top', bbox=props)

        plt.show()


def setup_interface(root):
    """
    Initializes Tk interface and widgets
    :return: NumericalIntegration object to use for exercises
    """
    # root.config(background='white')
    root.geometry("780x720")
    root.title('Function Approximation - Numerical Integration')

    # Give title to the window
    Label(root, text="Numerical Integration").pack()

    # Add a plot.
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(X_min, X_max)
    ax.set_ylim(Y_min, Y_max)
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()

    # Add the figure of matplotlib to the tkinter interface
    graph = FigureCanvasTkAgg(fig, master=root)
    graph.get_tk_widget().pack(side="top", fill='both')

    # Create object for numerical integration functions
    global numerical_integration
    numerical_integration = NumericalIntegration(root, graph, ax)

    # Prompt user for interval [a, b] and number of points using textbox and buttons.
    global a_text, b_text, n_text
    frame = Frame(root)
    frame.pack(padx=10, pady=10)

    a_text = Text(frame, width=20, height=2)
    a_text.grid(row=0, column=0, padx=10)

    button_a = Button(frame, text="Set value for 'a'", command=numerical_integration.set_a)
    button_a.grid(row=1, column=0, pady=10)

    b_text = Text(frame, width=20, height=2)
    b_text.grid(row=0, column=1, padx=10)

    button_b = Button(frame, text="Set value for 'b'", command=numerical_integration.set_b)
    button_b.grid(row=1, column=1, pady=10)

    n_text = Text(frame, width=20, height=2)
    n_text.grid(row=0, column=2)

    button_n = Button(frame, text="Set number of points", command=numerical_integration.set_nb_points)
    button_n.grid(row=1, column=2, pady=10)

    # Plot graph using initial data given.
    button_n = Button(frame, text="Plot function", command=numerical_integration.plot_function)
    button_n.grid(row=2, column=1, pady=10)

    # Buttons to start exercises
    button_n = Button(frame, text="Exercise 9", command=exercise9)
    button_n.grid(row=3, column=0, pady=10)
    button_n = Button(frame, text="Exercise 10", command=exercise10)
    button_n.grid(row=3, column=1, pady=10)

    # Button to display results
    button_res = Button(frame, text="Show results", command=numerical_integration.display_results)
    button_res.grid(row=3, column=2, pady=10)

    # Button to display results for last question
    button_res = Button(frame, text="Show last question results", command=numerical_integration.plot_last_question)
    button_res.grid(row=3, column=3, pady=10)


def reset_results():
    results['rectangle_inferieur'] = 0
    results['rectangle_superieur'] = 0
    results['trapeze'] = 0
    results['simpson'] = 0


def exercise9():
    reset_results()
    # Add and display exercise coordinates.

    # Rectangle inférieur method
    results['rectangle_inferieur'] = numerical_integration.rectangles_inferieurs("sin", 0, np.pi / 2, 4, None)

    # Rectangle supérieur method
    results['rectangle_superieur'] = numerical_integration.rectangles_superieurs("sin", 0, np.pi / 2, 4, None)

    # Simpson method
    results['simpson'] = numerical_integration.simpson("sin", 0, np.pi / 2, 4)

    # Rectangle trapèze method
    results['trapeze'] = numerical_integration.rectangles_trapeze("sin", 0, np.pi / 2, 4, None)


def exercise10():
    reset_results()

    y_coords = [30.0, 31.63, 33.44, 35.47, 37.75, 40.33, 43.29, 46.70, 50.67]

    # Simpson method
    results['simpson'] = numerical_integration.simpson_with_mapping(0, 80, 8, y_coords)

    # Rectangle trapèze method
    results['trapeze'] = numerical_integration.rectangles_trapeze(None, 0, 80, 8, y_coords)


def init():
    """
    Initialize the GUI for numerical integration
    :return:
    """
    root = Tk()
    setup_interface(root)

    # Test with functions from class exercises
    # exercise9()
    # exercise10()

    root.mainloop()


if __name__ == '__main__':
    init()
