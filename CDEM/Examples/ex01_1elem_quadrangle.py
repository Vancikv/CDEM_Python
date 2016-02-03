'''
Created on 8. 11. 2015

@author: Kapsak
'''

from CDEM import *
import numpy as np

# Force (n2x, n3x, n3y, n4x, n4y)
F_normal = np.array([0., 0., -0.5, 0., -0.5])
F_shear = np.array([0., 0.5, 0., 0.5, 0.])
F_mixed = F_shear + F_normal
tolerance = 10e-12
maxiter = 500
node_to_plot = None
time_step = .01

n1 = Node(x=0., y=0., supports=[1, 1])
n2 = Node(x=1., y=0., supports=[0, 1])
n3 = Node(x=1., y=1., F_ext=[0.5, 0.])
n4 = Node(x=0., y=1., F_ext=[0.5, 0.])

e1 = ElemQuadrangleLin(nodes=[1, 2, 3, 4], E=100., nu=0., density=1., thickness=1., alfaC=1.)

dom = DomainCDEMDR(c_n=1000., c_s=1000., elements=[e1], nodes=[n1, n2, n3, n4])

print '______NORMAL LOAD______'
n3.F_ext = F_normal[1:3]
n4.F_ext = F_normal[3:5]
e1.stiffdef = '2triangles'

dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2, 4, 5, 6, 7], :]
K = K[:, [2, 4, 5, 6, 7]]
Ki = np.linalg.inv(K)
print 'Displacements DR (STIFFNESS 2TRI):    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
print 'Displacements FEM (STIFFNESS 2TRI):    ', np.dot(Ki, F_normal)

e1.stiffdef = '4triangles'
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2, 4, 5, 6, 7], :]
K = K[:, [2, 4, 5, 6, 7]]
Ki = np.linalg.inv(K)
print 'Displacements DR (STIFFNESS 4TRI):    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
print 'Displacements FEM (STIFFNESS 4TRI):    ', np.dot(Ki, F_normal)

e1.stiffdef = 'isoparametric'
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2, 4, 5, 6, 7], :]
K = K[:, [2, 4, 5, 6, 7]]
Ki = np.linalg.inv(K)
print 'Displacements DR (STIFFNESS ISO):    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
print 'Displacements FEM (STIFFNESS ISO):    ', np.dot(Ki, F_normal), '\n\n'



print '______SHEAR LOAD______'
n3.F_ext = F_shear[1:3]
n4.F_ext = F_shear[3:5]
e1.stiffdef = '2triangles'
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2, 4, 5, 6, 7], :]
K = K[:, [2, 4, 5, 6, 7]]
Ki = np.linalg.inv(K)
print 'Displacements DR (STIFFNESS 2TRI):    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
print 'Displacements FEM (STIFFNESS 2TRI):    ', np.dot(Ki, F_shear)

e1.stiffdef = '4triangles'
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2, 4, 5, 6, 7], :]
K = K[:, [2, 4, 5, 6, 7]]
Ki = np.linalg.inv(K)
print 'Displacements DR (STIFFNESS 4TRI):    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
print 'Displacements FEM (STIFFNESS 4TRI):    ', np.dot(Ki, F_shear)

e1.stiffdef = 'isoparametric'
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2, 4, 5, 6, 7], :]
K = K[:, [2, 4, 5, 6, 7]]
Ki = np.linalg.inv(K)
print 'Displacements DR (STIFFNESS ISO):    ', np.hstack((n2.v_disp[0], n3.v_disp, n4.v_disp))
print 'Displacements FEM (STIFFNESS ISO):    ', np.dot(Ki, F_shear)
