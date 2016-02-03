'''
Created on 8. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *
'''
n1 = Node(x=0., y=0., supports=[1, 1])
n2 = Node(x=1., y=0., supports=[0, 1])
n3 = Node(x=1., y=1., neighbors=[0, 6])
n4 = Node(x=0., y=1., neighbors=[5, 0])

n5 = Node(x=0., y=1., neighbors=[0, 4])
n6 = Node(x=1., y=1., neighbors=[3, 0])
n7 = Node(x=1., y=2., F_ext=[2.5, 0.])
n8 = Node(x=0., y=2., F_ext=[2.5, 0.])

e1 = ElemQuadrangleLin(nodes=[1, 2, 3, 4], E=100., nu=0., density=1., thickness=1., alfaC=15.)
e2 = ElemQuadrangleLin(nodes=[5, 6, 7, 8], E=100., nu=0., density=1., thickness=1., alfaC=15.)

dom = DomainCDEMStatic(c_n=500., c_s=100., elements=[e1, e2], nodes=[n1, n2, n3, n4, n5, n6, n7, n8])

dom.solve()
dom.plot(magnitude=1.0)
'''
dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             c_n=100000., c_s=100000., domtype='CDEMStatic',
                             load=[((35., 65.), (25., 35.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

dom.solve()
dom.plot(magnitude=50.)

