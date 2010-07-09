#! /usr/bin/env python

from utils import conv_graph2
from hydrogen_romanowski import R_x as R2_x, R_y as R2_y
from hydrogen_certik_rR import R_x, R_y
from hydrogen_certik_R import R_x as R3_x, R_y as R3_y

conv_graph2(R_x, R_y, R2_x, R2_y, 0)
conv_graph2(R_x, R_y, R2_x, R2_y, 1)
conv_graph2(R_x, R_y, R2_x, R2_y, 2)
