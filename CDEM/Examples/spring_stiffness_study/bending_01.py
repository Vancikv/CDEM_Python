'''
Created on 10. 11. 2015

@author: Kapsak

A beam subjected to 3point bending.
FEM solution is compared with static CDEM solution with contact stiffness
varying by powers of 10.
Default contact stiffness for static CDEM solution is obtained
as the largest number across all the local stiffness matrices.
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np
import tabulate as tbl
from CDEM.file_handler import get_outfile
headers = ['Stiffness factor order', 'c_n', 'c_s',
                            'w_u', 'w_u/w_FEM', 's_u', 's_u/s_FEM',
                            'w_s', 'w_s/w_FEM', 's_s', 's_s/s_FEM',
                            'w_FEM', 's_FEM', 'Longest edge', 'Shortest edge']
table = {key:[] for key in headers}
table['Stiffness factor order'] = range(-1, 8)
dom1 = get_2D_quadrangle_FEM_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             load=[((45., 55.), (25., 35.), [0., -100.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

dom1.solve()
# dom1.plot(magnitude=50.)
w_FEM = dom1.nodes[23].v_disp[1]
s_FEM = dom1.get_stress_norm([0])
print 'Reference node coordinates: x = %f, y = %f' % (dom1.nodes[23].x, dom1.nodes[23].y)
print 'FEM deflection: %.6f' % w_FEM
print 'FEM sigma_x L2 norm: %.6f' % s_FEM
dom1.plot(magnitude=50., savefile=get_outfile(['Dropbox', 'CDEM', 'spring_stiffness_study'], 'beam01.png'), show=False)

dom2 = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=1., domtype='CDEMStatic',
                             load=[((45., 55.), (25., 35.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

C_def = dom2.C
dom2.solve()
minl = min([min(nd.l_edges) for nd in dom2.nodes])
maxl = max([max(nd.l_edges) for nd in dom2.nodes])
f = open(get_outfile(['Dropbox', 'CDEM', 'spring_stiffness_study'], 'beam01_stiffness_variation.txt'), 'w')
for i in table['Stiffness factor order']:
    dom2.C = C_def * (10 ** i)
    dom2.stiffness_scaling = False
    table['w_FEM'].append(w_FEM)
    table['s_FEM'].append(s_FEM)
    table['Longest edge'].append(maxl)
    table['Shortest edge'].append(minl)
    dom2.solve()
    # dom2.plot(magnitude=50., savefile='c:\\Users\\Werner\\Dropbox\\CDEM\\spring_stiffness_study\\bending01\\stiff_fact_%d.png' % i, show=False)
    table['c_n'].append(dom2.C[0, 0])
    table['c_s'].append(dom2.C[1, 1])
    # f.write('Reference node coordinates: x = %f, y = %f\n' % (dom2.nodes[58].x, dom2.nodes[58].y))
    table['w_u'].append(dom2.nodes[58].v_disp[1])
    table['w_u/w_FEM'].append(dom2.nodes[58].v_disp[1] / w_FEM)
    table['s_u'].append(dom2.get_stress_norm([0]))
    table['s_u/s_FEM'].append(dom2.get_stress_norm([0]) / s_FEM)
    dom2.stiffness_scaling = True
    dom2.solve()
    table['w_s'].append(dom2.nodes[58].v_disp[1])
    table['w_s/w_FEM'].append(dom2.nodes[58].v_disp[1] / w_FEM)
    table['s_s'].append(dom2.get_stress_norm([0]))
    table['s_s/s_FEM'].append(dom2.get_stress_norm([0]) / s_FEM)
texttable = tbl.tabulate(zip(*[table[col] for col in headers]), headers=headers, tablefmt="orgtbl")
f.write(texttable)
f.close()

