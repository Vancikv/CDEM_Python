'''
Created on 10. 11. 2015

@author: Kapsak

'''

from CDEM import *
import numpy as np

dom = get_2D_quadrangle_domain(ni=8, nj=96, lx=2., ly=0.2, E=25000000000., nu=0.2, density=2500.,
                             thickness=0.2, stiffdef='isoparametric',
                             spring_stiffness_factor=10.0, domtype='CDEMStatic',
                             load=[((0.99999, 1.00001), (0.19999, 0.20001), [0., -200000.])],
                             supports=[((-0.00001, 0.00001), (-0.00001, 0.00001), [1, 1]), ((1.99999, 2.00001), (-0.00001, 0.00001), [0, 1])])

dom.save_to_text_file("c:/Users/Kapsak/Dropbox/CDEM/Cpp_data/3pb_binary/3pb.txt")
dom.solve()
dom.plot(magnitude=10.)
