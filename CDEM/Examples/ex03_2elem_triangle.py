'''
Created on 8. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *

n1 = Node(x=0., y=0., supports=[1, 1], neighbors=[6, 0])
n2 = Node(x=1., y=0., supports=[0, 1])
n3 = Node(x=1., y=1., F_ext=[0., -0.5], neighbors=[0, 4])
n4 = Node(x=1., y=1., neighbors=[3, 0])
n5 = Node(x=0., y=1., F_ext=[0., -0.5])
n6 = Node(x=0., y=0., neighbors=[0, 1])

e1 = ElemTriangleLin(nodes=[1, 2, 3], E=100., nu=0., density=1., thickness=1., alfaC=20.)
e2 = ElemTriangleLin(nodes=[4, 5, 6], E=100., nu=0., density=1., thickness=1., alfaC=20.)

dom = Domain(c_n=500., c_s=500., elements=[e1, e2], nodes=[n1, n2, n3, n4, n5, n6])
dom.solve(dt=0.005, tol=0.00000001, maxiter=200, node_to_plot=3)
dom.plot(magnitude=5.)

print 'DR displacements from shear load (n2x, n3x, n3y, n5x, n5y):\n', np.hstack((n2.v_disp[0], n3.v_disp, n5.v_disp))
