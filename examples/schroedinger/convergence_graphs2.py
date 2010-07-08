from math import e, log, exp

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel


Rn0_x = [25, 31, 37, 43, 49, 67, 79, 103, 115, 133]
#Rn1_x = [23, 35, 41, 47, 59, 71, 83, 95, 107]
#Rn2_x = [23, 29, 35, 41, 47, 65, 77, 95, 107]

R10_y = [0.16676173689459373, 0.038254689188603141, 0.0010677432574189294,
        2.2931902048939357e-06, 1.0921685889009325e-08, 9.8654451274882149e-10,
        9.3617336105467075e-12, 7.0984329525458634e-12, 6.4437344349244086e-13,
        2.8999025403209089e-13]

R20_y = [0.028966792038071301, 0.0030262717329101185, 1.6489249960768837e-05,
        1.7408326863366241e-08, 9.1130985258036645e-09, 3.8598236690168264e-09,
        1.5295875677168169e-11, 2.8138602559124593e-13, 1.0014211682118912e-13,
        -2.3553381467422696e-13]

R30_y = [0.010422065707353885, 0.00068778277817766603, 1.7278464779918457e-06,
        2.8916446515037109e-08, 2.843261607404779e-08, 8.0352775128078591e-09,
        2.276666369316338e-09, 5.9793628381932251e-12, -1.602190602412179e-14,
        4.0821512836686225e-14]

R40_y = [0.0050153980309137862, 0.00025054408623355687, 4.5227925233801147e-07,
        5.223349718377901e-08, 5.2148303317928191e-08, 1.2370829003544026e-08,
        7.4851823489430203e-09, 5.0776650245554222e-11, 5.7518920182353384e-13,
        5.4471357979757329e-13]

R21_y = [0.018109721060183215, 0.0044070390445993674, 7.1159741065113247e-05,
        8.2037840076765178e-08, 3.2771956398613966e-09, 9.2940505291849718e-10,
        1.1640521879741073e-11, 1.5290546606649968e-12, 7.2601646916581331e-13,
        1.0146050666293149e-13, 6.8389738316909643e-14,
        -2.2787327580431338e-14, -7.1609385088322597e-15,
        8.7735374521002996e-14, -3.4555691641457997e-14]
R31_y = [0.0099474675951749178, 0.0018084849047237742, 1.4646370057380675e-05,
        1.2492772066829971e-08, 4.2047221121555012e-09, 5.9871994861904199e-10,
        2.0967942659932959e-11, 1.9885343371939257e-12, 5.8615612363865921e-13,
        2.7159524629283283e-13, 3.1863400806741993e-14,
        -6.1825544683813405e-15, 1.366962099069724e-15,
        -1.5742268599794329e-13, -4.3069714461552167e-14]
R41_y = [0.0057412836805739824, 0.00088112116982369884, 4.8025338958666841e-06,
        1.9832536838487735e-08, 1.7926716239929918e-08, 5.8603496569387126e-09,
        5.550401666054583e-09, 7.9150505594149934e-12, 8.2727227845857954e-13,
        1.0153683449587447e-12, 5.2366791458702266e-13,
        -4.3728909382423353e-14, -7.8770323597154857e-14,
        1.4579309981499478e-13, 1.9487536584428256e-13]
R51_y = [0.0035449387734955082, 0.00049043403644849146, 2.1651823179751062e-06,
        8.1484721182328856e-08, 8.0809502499279429e-08, 5.1479379131053049e-08,
        5.1303198862762134e-08, 8.5330239402159336e-11, 2.9134542001152397e-12,
        2.5417064286603619e-12, 2.6824237275846485e-12, 1.0858570986815863e-12,
        7.1252032052271375e-14, -2.5522639557351567e-13,
        -2.3429175266542757e-14]

#R32_y
#R42_y
#R52_y
#R62_y

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

figure()
do_plot(Rn0_x, R10_y, 1, 0)
do_plot(Rn0_x, R20_y, 2, 0)
do_plot(Rn0_x, R30_y, 3, 0)
do_plot(Rn0_x, R40_y, 4, 0)
savefig("conv_dof_l_0.png")

#figure()
#do_plot(Rn1_x, R21_y, 2, 1)
#do_plot(Rn1_x, R31_y, 3, 1)
#do_plot(Rn1_x, R41_y, 4, 1)
#do_plot(Rn1_x, R51_y, 5, 1)
#savefig("conv_dof_l_1.png")
#
#figure()
#do_plot(Rn2_x, R32_y, 3, 2)
#do_plot(Rn2_x, R42_y, 4, 2)
#do_plot(Rn2_x, R52_y, 5, 2)
#do_plot(Rn2_x, R62_y, 6, 2)
#savefig("conv_dof_l_2.png")
#
#show()
