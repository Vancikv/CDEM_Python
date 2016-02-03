'''
Created on Nov 27, 2015

@author: Werner
'''

from CDEM.basic import *

# Force (n2x, n3x, n3y, n4x, n4y)
tolerance = 1e-8
maxiter = 100
maxiniter = 1000
node_to_plot = None
time_step = .01

n1 = Node(x=0., y=0., supports=[1, 1])
n2 = Node(x=1., y=0., supports=[0, 1])
n3 = Node(x=1., y=1., F_ext=[ 0., -0.5])
n4 = Node(x=0., y=1., F_ext=[ 0., -0.5])

e1 = ElemQuadrangleLin(nodes=[1, 2, 3, 4], E=100., nu=0., density=1., thickness=1., alfaC=1.)

dom = DomainCDEMDR(c_n=1000., c_s=1000., elements=[e1], nodes=[n1, n2, n3, n4])

dom.solve_kin_damping(dt=time_step, tol=tolerance, maxiter=maxiter, maxiniter=maxiniter, verbose=False)
print 'Displacements kinematic damping:    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
