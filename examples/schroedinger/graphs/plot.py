#! /usr/bin/env python

from utils import conv_graph2
import hydrogen_romanowski as romanowski
import hydrogen_pask_R as pask
import hydrogen_certik_rR
import hydrogen_certik_R

conv_graph2(pask.R_x, pask.R_y, romanowski.R_x, romanowski.R_y, 0, "R", eigs=3)
#conv_graph2(R_x, R_y, R2_x, R2_y, 0, "rR")
#conv_graph2(R_x, R_y, R2_x, R2_y, 1, "rR")
#conv_graph2(R_x, R_y, R2_x, R2_y, 2, "rR")

#conv_graph2(R3_x, R3_y, R_x, R_y, 0, "R")
#conv_graph2(R3_x, R3_y, R_x, R_y, 1, "R")
#conv_graph2(R3_x, R3_y, R_x, R_y, 2, "R")
