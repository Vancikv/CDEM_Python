'''
Created on Jan 6, 2016

@author: vancik
'''

from CDEM.file_handler import get_outfile
from CDEM import *

# dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
#                              thickness=1., stiffdef='isoparametric',
#                              spring_stiffness_factor=10., domtype='CDEMStatic',
#                              load=[((49, 51.), (29., 31.), [0., -50.])],
#                              supports=[((-1., 1.), (-1., 1.), [1, 1]), ((99., 101.), (-1., 1.), [0, 1])])

dom = get_2D_quadrangle_domain(ni=3, nj=20, lx=2., ly=0.2, E=25000000000., nu=0.2, density=2500., alfaC=0.0,
                             thickness=0.2, stiffdef='isoparametric',
                             spring_stiffness_factor=10., domtype='CDEMDR',
                             load=[((0.99, 1.01), (0.19, 0.21), [0., -200000.])],
                             supports=[((-0.01, 0.01), (-0.01, 0.01), [1, 1]), ((1.99, 2.01), (-0.01, 0.01), [0, 1])])

f = open(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input'], '3point_bending_real.in'), 'w')

for el in dom.elements:
    el.calc_normal_vectors()
springs = []
nvects = []
for i, node in enumerate(dom.nodes):
    for nb, nv in zip(node.neighbors, node.n_vects):
        if nb != 0:
            couple = (min(i + 1, nb), max(i + 1, nb))
            if couple not in springs:
                springs.append(couple)
                nvects.append(nv)

f.write('ndofman %d nelem %d\n' % (len(dom.nodes), len(dom.elements) + 2 * len(springs)))
for i, node in enumerate(dom.nodes):
    f.write('Node %d coords 2 %g %g\n' % (i + 1, node.x, node.y))
for i, el in enumerate(dom.elements):
    nds = el.nodes
    f.write('PlaneStress2d %d nodes 4 %d %d %d %d\n' % (i + 1, nds[0], nds[1], nds[2], nds[3]))
for i, sp in enumerate(springs):
    nv = nvects[i]
    f.write('Spring %d nodes 2 %d %d mode 1 k %f orientation 3 %f %f %f\n' % (len(dom.elements) + 1 + 2 * i, sp[0], sp[1], dom.C[0, 0], nv[0], nv[1], 0.))
    f.write('Spring %d nodes 2 %d %d mode 1 k %f orientation 3 %f %f %f\n' % (len(dom.elements) + 1 + 2 * i + 1, sp[0], sp[1], dom.C[1, 1], nv[1], -nv[0], 0.))

f.close()
