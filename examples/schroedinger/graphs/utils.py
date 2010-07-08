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
