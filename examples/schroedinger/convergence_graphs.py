from math import e, log, exp

from pylab import plot, show, savefig, grid, gca, legend

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



g = Convert([
    (274, 1),
    (845, 1e-12),
], scale="log", axis="y")

Rn0_x = [23, 29, 41, 47, 53, 65, 83, 107]
R10_y = g.convert([
    (254, 305),
    (319, 325),
    (382, 376),
    (447, 474),
    (512, 611),
    (575, 723),
    (640, 768),
    (703, 822),
    ])
R20_y = g.convert([
    (254, 338),
    (319, 367),
    (383, 443),
    (447, 569),
    (511, 679),
    (575, 715),
    (639, 832),
    ])
R30_y = g.convert([
    (255, 358),
    (319, 393),
    (383, 480),
    (448, 616),
    (512, 681),
    (576, 707),
    (640, 835),
    ])

R40_y = g.convert([
    (254, 372),
    (320, 410),
    (383, 504),
    (448, 635),
    (512, 653),
    (575, 655),
    (640, 787),
    ])

def do_plot(x, y, label, marker="o"):
    z = zip(x, y)
    x = [_[0] for _ in z]
    y = [_[1] for _ in z]
    plot(x, y, "k-")
    plot(x, y, "k" + marker, label=label)

do_plot(Rn0_x, R10_y, "$R_{10}$", "s")
do_plot(Rn0_x, R20_y, "$R_{20}$", "o")
do_plot(Rn0_x, R30_y, "$R_{30}$", "^")
do_plot(Rn0_x, R40_y, "$R_{40}$", "v")



grid(True)
ax = gca()
ax.set_yscale("log")
legend()
show()
#savefig("dofs_l_2.png")
