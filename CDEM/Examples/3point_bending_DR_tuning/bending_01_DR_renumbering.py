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

headers = ['d_t', 'n_iter', 'density', 'w_end', 's_end', 'C_min', 'C_max', 'alfaC', 'stiffness_factor']
table = {key:[] for key in headers}
table['d_t'] = 2 * [0.001]
table['density'] = 2 * [1.]
table['alfaC'] = 2 * [50 * ((25.*1.) ** 0.5)]
table['stiffness_factor'] = 2 * [100.]
n_iter = 1

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1., alfaC=1.0,
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=1., domtype='CDEMDR',
                             load=[((45., 55.), (25., 35.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

C_def = dom.C
dom.solve(dt=0.01, tol=1e-16, maxiter=1)

ro = table['density'][0]
alfaC = table['alfaC'][0]
sfact = table['stiffness_factor'][0]
dt = table['d_t'][0]
for el in dom.elements:
    el.density = ro
    el.alfaC = alfaC
for nd in dom.nodes:
    nd.v_velo = np.zeros(2)
    nd.v_disp = np.zeros(2)

dom.C *= sfact
dom.solve(dt=dt, tol=1e-16, maxiter=n_iter, node_to_plot=59, savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_tuning'], 'graph_ro=%.1f_dt=%.3f_alfaC=%f_stifffact=%.1f.png' % (ro, dt, alfaC, sfact)), show=False)
C_min = min([min([el.C[i, i] for i in range(el.C.shape[0])]) for el in dom.elements])
C_max = max([np.max(el.C) for el in dom.elements])
dom.plot(magnitude=100., savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_tuning'], 'deformed_ro=%.1f_dt=%.3f_alfaC=%f_stifffact=%.1f.png' % (ro, dt, alfaC, sfact)), show=False)
table['w_end'].append(dom.nodes[58].v_disp[1])
table['s_end'].append(dom.get_stress_norm([0]))
table['C_min'].append(C_min)
table['C_max'].append(C_max)
table['n_iter'].append(n_iter)
dom.C /= sfact

dom.renumber_elements()
for el in dom.elements:
    el.density = ro
    el.alfaC = alfaC
for nd in dom.nodes:
    nd.v_velo = np.zeros(2)
    nd.v_disp = np.zeros(2)

dom.C *= sfact
dom.solve(dt=dt, tol=1e-16, maxiter=n_iter, node_to_plot=59, savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_renumbering'], 'graph_renumbered_ro=%.1f_dt=%.3f_alfaC=%f_stifffact=%.1f.png' % (ro, dt, alfaC, sfact)), show=False)
C_min = min([min([el.C[i, i] for i in range(el.C.shape[0])]) for el in dom.elements])
C_max = max([np.max(el.C) for el in dom.elements])
dom.plot(magnitude=100., savefile=get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_renumbering'], 'deformed_renumbered_ro=%.1f_dt=%.3f_alfaC=%f_stifffact=%.1f.png' % (ro, dt, alfaC, sfact)), show=False)
table['w_end'].append(dom.nodes[58].v_disp[1])
table['s_end'].append(dom.get_stress_norm([0]))
table['C_min'].append(C_min)
table['C_max'].append(C_max)
table['n_iter'].append(n_iter)
dom.C /= sfact

f = open(get_outfile(['Dropbox', 'CDEM', '3point_bending_DR_renumbering'], '3point_bending_DR_renumbering_study.txt'), 'w')
texttable = tbl.tabulate(zip(*[table[col] for col in headers]), headers=headers, tablefmt="orgtbl")
f.write(texttable)
f.close()
