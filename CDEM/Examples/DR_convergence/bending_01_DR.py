'''
Created on 10. 11. 2015

@author: Kapsak

A beam subjected to 3point bending.
The contact stiffness factor is chosen as 100 in order
to achieve sufficient accuracy with respect to the FEM
solution.
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np
import tabulate as tbl
from CDEM.file_handler import get_outfile

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1., alfaC=1.0,
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=1., domtype='CDEMDR',
                             load=[((45., 55.), (25., 35.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])


C_def = dom.C
dom.solve(dt=0.01, tol=1e-16, maxiter=1)

headers = ['d_t', 'n_iter', 'density', 'w_end', 's_end', 'C_min', 'C_max', 'alfaC', 'stiffness_factor']
table = {key:[] for key in headers}
max_omega = (np.max(C_def)*10./25.) ** 0.5
recommended_step = 2. / max_omega
table['d_t'] = 4*[0.001]
table['density'] = 4 * [1.]
table['alfaC'] = 4 * [0.]
table['stiffness_factor'] = [10., 10., 10., 10.]
table['n_iter'] = [1000,2000,4000,6000]

for dt, ro, alfaC, sfact, n_iter in zip(table['d_t'], table['density'], table['alfaC'], table['stiffness_factor'], table['n_iter']):
    for el in dom.elements:
        el.density = ro
        el.alfaC = alfaC
    for nd in dom.nodes:
        nd.v_velo = np.zeros(2)
        nd.v_disp = np.zeros(2)
    dom.C *= sfact
    dom.solve(dt=dt, tol=1e-16, maxiter=n_iter, node_to_plot=59, savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_convergence'], 'graph_niter=%d.png' % n_iter), show=False)
    C_min = min([min([el.C[i, i] for i in range(el.C.shape[0])]) for el in dom.elements])
    C_max = max([np.max(el.C) for el in dom.elements])
    dom.plot(magnitude=100., savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_convergence'], 'deformed_niter=%d.png' % n_iter), show=False)
    table['w_end'].append(dom.nodes[58].v_disp[1])
    table['s_end'].append(dom.get_stress_norm([0]))
    table['C_min'].append(C_min)
    table['C_max'].append(C_max)
    table['n_iter'].append(n_iter)
    dom.C /= sfact
f = open(get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_convergence'], '3point_bending_DR_convergence.txt'), 'w')
texttable = tbl.tabulate(zip(*[table[col] for col in headers]), headers=headers, tablefmt="orgtbl")
f.write(texttable)
f.close()
