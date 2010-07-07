from pylab import plot, show, savefig, grid, gca

R52 = [5*1e-3, 1e-3, 1e-4, 1e-6, 3*1e-8, 1e-9, 7*1e-11, 5*1e-12, 3*1e-12]
R62 = [5*1e-3, 1e-3, 1e-4, 1e-6, 5*1e-8, 1e-8, 1e-10, 5*1e-12, 3*1e-12]

steps = 9

x = []
y = []
for n in range(steps):
    x.append(n+1)
    y.append(R52[n])
plot(x, y)


grid(True)
ax = gca()
ax.set_yscale("log")
savefig("dofs_l_2.png")
