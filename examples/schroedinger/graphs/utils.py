from math import log

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

def do_plot(x, y, n, l, color="k"):
    n_r = n - l - 1
    styles = {0: "-s", 1: "--o", 2: ":^", 3: "-.v"}
    plot(x, y, color + styles[n_r], label="$R_{%d%d}$" % (n, l))

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$")
    ax.set_yscale("log")
    title("l=%d" % l)
    legend()

def conv_graph(R_x, R_y, l):
    filename = "conv_dof_l_%d.png" % l
    print "Creating %s" % filename
    figure()
    for n in range(l+1, l+5):
        do_plot(R_x[l], R_y[(n, l)], n, l)
    savefig(filename)

class Convert(object):

    def __init__(self, points, scale="linear", axis="x"):
        if scale == "log":
            self._x1, self._val1 = points[0]
            self._x2, self._val2 = points[1]
        else:
            raise NotImplementedError("Scale not implemented yet")

        if axis == "x":
            self._axis = 0
        elif axis == "y":
            self._axis = 1
        else:
            raise ValueError("axis must be either 'x' or 'y'")

    def conv(self, x):
        # transform 'x' to the interval [0, 1]:
        p01 = float((x - self._x1))/(self._x2 - self._x1)
        # transform 'x' to the interval [v1, v2]:
        v1 = log(self._val1)/log(10)
        v2 = log(self._val2)/log(10)
        pv1v2 = p01 * (v2 - v1) + v1
        # finally just exponentiate it:
        return 10**pv1v2

    def convert(self, data):
        return [self.conv(x[self._axis]) for x in data]


