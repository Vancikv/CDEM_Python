'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np
import tabulate as tbl
from CDEM.file_handler import get_outfile

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1., alfaC=50.*(25.**0.5),
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=100., domtype='CDEMStatic',
                             load=[((45., 55.), (25., 35.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])


dom.solve(timelog=get_outfile(['Dropbox', 'CDEM', '3point_bending_timing'], '3point_bending_Static_timing.txt'))
