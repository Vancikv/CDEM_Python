'''
Created on 10. 11. 2015

@author: Kapsak

'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np
import tabulate as tbl
from CDEM.file_handler import get_outfile

dom = get_2D_quadrangle_domain(ni=3, nj=20, lx=2., ly=0.2, E=25000000000., nu=0.2, density=2500., alfaC=0.0,
                             thickness=0.2, stiffdef='isoparametric',
                             spring_stiffness_factor=10., domtype='CDEMDR',
                             load=[((0.99, 1.01), (0.19, 0.21), [0., -200000.])],
                             supports=[((-0.01, 0.01), (-0.01, 0.01), [1, 1]), ((1.99, 2.01), (-0.01, 0.01), [0, 1])])


headers = ['d_t', 'n_iter', 'density', 'w_end', 's_end', 'C_min', 'C_max', 'alfaC', 'stiffness_factor']
table = {key:[] for key in headers}

table['n_iter'] = [20000, 40000]

for n_iter in table['n_iter']:
    for nd in dom.nodes:
        nd.v_velo = np.zeros(2)
        nd.v_disp = np.zeros(2)
    dom.solve(t_load=0.004, t_max=0.004, tol=1e-16, maxiter=n_iter, node_to_plot=119, savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'graph_niter=%.1f.png' % n_iter),
              savetextfile=get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=%.1f.txt' % (n_iter)), show=False,
              output_freq=1000)
    dom.plot(magnitude=100., savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'deformed_niter=%.1f.png' % (n_iter)), show=False)
