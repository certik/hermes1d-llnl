from math import e, log, exp

from pylab import plot, show, savefig, grid, gca, legend, figure

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
Rn1_x = [23, 29, 35, 41, 47, 65, 77, 95, 107]

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

R32_y = g.convert([
    (973, 274),
    (1037, 339),
    (1102, 377),
    (1165, 479),
    (1230, 635),
    (1295, 711),
    (1359, 807),
    (1425, 809),
    (1488, 828),
])

R42_y = g.convert([
    (973, 278),
    (1038, 335),
    (1102, 384),
    (1166, 493),
    (1230, 619),
    (1294, 698),
    (1359, 785),
    (1423, 812),
    (1488, 818),
])

R52_y = g.convert([
    (972, 288),
    (1038, 341),
    (1103, 391),
    (1166, 506),
    (1231, 602),
    (1295, 669),
    (1359, 752),
    (1424, 818),
    (1488, 825),
])

R62_y = g.convert([
    (973, 296),
    (1037, 346),
    (1102, 401),
    (1166, 516),
    (1231, 581),
    (1295, 600),
    (1360, 725),
    (1424, 811),
    (1488, 829),
])

def do_plot(x, y, label, marker="o"):
    z = zip(x, y)
    x = [_[0] for _ in z]
    y = [_[1] for _ in z]
    plot(x, y, "k-")
    plot(x, y, "k" + marker, label=label)

    grid(True)
    ax = gca()
    ax.set_yscale("log")
    legend()

do_plot(Rn0_x, R10_y, "$R_{10}$", "s")
do_plot(Rn0_x, R20_y, "$R_{20}$", "o")
do_plot(Rn0_x, R30_y, "$R_{30}$", "^")
do_plot(Rn0_x, R40_y, "$R_{40}$", "v")

figure()

do_plot(Rn1_x, R32_y, "$R_{32}$", "s")
do_plot(Rn1_x, R42_y, "$R_{42}$", "o")
do_plot(Rn1_x, R52_y, "$R_{52}$", "^")
do_plot(Rn1_x, R62_y, "$R_{62}$", "v")



show()
#savefig("dofs_l_2.png")
