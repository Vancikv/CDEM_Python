'''
Created on 10. 11. 2015

@author: Kapsak

Beam from bending04 with a different mesh.
'''

from CDEM.basic import *
from CDEM.utils import *
from CDEM.file_handler import get_outfile
import numpy as np
import tabulate as tbl

headers = ['Stiffness factor order', 'c_n', 'c_s',
                            'w_u', 'w_u/w_FEM', 's_u', 's_u/s_FEM',
                            'w_s', 'w_s/w_FEM', 's_s', 's_s/s_FEM',
                            'w_FEM', 's_FEM', 'Longest edge', 'Shortest edge']
table = {key:[] for key in headers}
table['Stiffness factor order'] = range(-1, 8)


dom2 = get_uneven_mesh(ni0=3, nj=16, lx=100., ly=30., stencil_arr=['0', '0', '0', 'f', '0', 'f', '0', '0', '0', '0', 'c', '0', 'c', '0', '0', '0'],
                             E=10000., nu=0.2, density=1.,
                             thickness=1., alfaC=1., stiffdef='isoparametric',
                             spring_stiffness_factor=1., domtype='CDEMStatic',
                             load=[((45.0, 55.0), (29.9, 30.1), [0., -50. / 8.])],
                             supports=[((-0.1, 0.1), (-0.1, 0.1), [1, 1]), ((99.9, 100.1), (-0.1, 0.1), [0, 1])])

dom1 = convert_CDEM_to_FEM_domain(dom2)
for nd in dom1.nodes:
    nd.F_ext[0] *= 200. / 87.5
    nd.F_ext[1] *= 200. / 87.5

dom1.solve()
# dom1.plot(magnitude=50.)
refnodes = dom1.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
w_FEM = dom1.nodes[refnd].v_disp[1]
s_FEM = dom1.get_stress_norm([0])
print 'Reference node coordinates: x = %f, y = %f' % (dom1.nodes[refnd].x, dom1.nodes[refnd].y)
print 'FEM deflection: %.6f' % w_FEM
print 'FEM sigma_x L2 norm: %.6f' % s_FEM
dom1.plot(magnitude=50., savefile=get_outfile(['Dropbox', 'CDEM', 'spring_stiffness_study'], 'beam05.png'), show=False)

C_def = dom2.C
dom2.solve()
minl = min([min(nd.l_edges) for nd in dom2.nodes])
maxl = max([max(nd.l_edges) for nd in dom2.nodes])
refnodes = dom2.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]

from CDEM.file_handler import get_outfile
f = open(get_outfile(['Dropbox', 'CDEM', 'spring_stiffness_study'], 'beam05_stiffness_variation.txt'), 'w')
for i in table['Stiffness factor order']:
    dom2.C = C_def * (10 ** i)
    dom2.stiffness_scaling = False
    table['w_FEM'].append(w_FEM)
    table['s_FEM'].append(s_FEM)
    table['Longest edge'].append(maxl)
    table['Shortest edge'].append(minl)
    dom2.solve()
    table['c_n'].append(dom2.C[0, 0])
    table['c_s'].append(dom2.C[1, 1])
    table['w_u'].append(dom2.nodes[refnd].v_disp[1])
    table['w_u/w_FEM'].append(dom2.nodes[refnd].v_disp[1] / w_FEM)
    table['s_u'].append(dom2.get_stress_norm([0]))
    table['s_u/s_FEM'].append(dom2.get_stress_norm([0]) / s_FEM)
    dom2.stiffness_scaling = True
    dom2.solve()
    table['w_s'].append(dom2.nodes[refnd].v_disp[1])
    table['w_s/w_FEM'].append(dom2.nodes[refnd].v_disp[1] / w_FEM)
    table['s_s'].append(dom2.get_stress_norm([0]))
    table['s_s/s_FEM'].append(dom2.get_stress_norm([0]) / s_FEM)
texttable = tbl.tabulate(zip(*[table[col] for col in headers]), headers=headers, tablefmt="orgtbl")
f.write(texttable)
f.close()
print texttable

