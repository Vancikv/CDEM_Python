'''
Created on 24. 11. 2015

@author: Kapsak
'''
from CDEM.basic import *

n1 = Node(x=0., y=0., supports=[1, 1])
n2 = Node(x=1., y=0., supports=[0, 1])
n3 = Node(x=1., y=1., F_ext=[0.5, 0.])
n4 = Node(x=0., y=1., F_ext=[0.5, 0.])

e1 = ElemQuadrangleLin(nodes=[1, 2, 3, 4], E=100., nu=0., density=1., thickness=1., alfaC=15.)

dom = DomainCDEMStatic(c_n=1000., c_s=1000., elements=[e1], nodes=[n1, n2, n3, n4])

e1.set_matrices()
print e1.K
n2.x = 10.
n3.x = 10.
e1.set_matrices()
print e1.K

n3.y = 10.
n4.y = 10.
e1.set_matrices()
print e1.K
