'''
Created on Nov 24, 2015

@author: vancik
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np

dom1 = get_uneven_mesh(ni0=2, nj=2, lx=30., ly=30., stencil_arr=['f', 'f'],
                             E=10000., nu=0.2, density=1.,
                             thickness=1., alfaC=1., stiffdef='isoparametric',
                             c_n=100000., c_s=100000.,
                             spring_stiffness_factor=1., domtype='CDEMStatic',
                             load=[((49.9, 50.1), (29.9, 30.1), [0., -50.])],
                             supports=[((-0.1, 0.1), (-0.1, 0.1), [1, 1]), ((99.9, 100.1), (-0.1, 0.1), [0, 1])])

# dom1.solve()
dom1.plot(magnitude=5.)
els = dom1.elements
nds = dom1.nodes

'''
for i, el in enumerate(els):
    print 'Element %d:' % (i + 1)
    for n in el.nodes:
        print '    Node %d ... neighbors [%d, %d]' % (n, nds[n - 1].neighbors[0], nds[n - 1].neighbors[1])


isclose = lambda n1, n2, tol = 1e-6: (((n2[0] - n1[0]) ** 2 + (n2[1] - n1[1]) ** 2) ** 0.5) < tol

nds = dom1.nodes
for nd in nds:
    #print nd.neighbors
    for n in nd.neighbors:
        if n != 0:
            if not isclose((nds[n-1].x,nds[n-1].y),(nd.x,nd.y)):
                print 'node:', (nd.x,nd.y)
                print 'wrong neighbor:', (nds[n-1].x,nds[n-1].y)
'''
