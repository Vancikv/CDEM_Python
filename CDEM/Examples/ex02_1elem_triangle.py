'''
Created on 8. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *

F_normal = np.array([0., 0., -0.5])
F_shear = np.array([0., 0.5, 0.])
F_mixed = F_shear + F_normal
tolerance = 10e-8
maxiter = 200
node_to_plot = None
time_step = .01

n1 = Node(x=0., y=0., supports=[1, 1])
n2 = Node(x=1., y=0., supports=[0, 1])
n3 = Node(x=1., y=1., F_ext=F_normal[1:])

e1 = ElemTriangleLin(nodes=[1, 2, 3], E=100., nu=0., density=1., thickness=1., alfaC=15.)

dom = Domain(c_n=1000., c_s=1000., elements=[e1], nodes=[n1, n2, n3])
print '-------- nu = 0.0 --------\n\n'
# Normal load
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2,4,5],:]
K = K[:,[2,4,5]]
Ki = np.linalg.inv(K)
print 'FEM displacements from normal load (n2x, n3x, n3y):\n', np.dot(Ki, F_normal)
print 'DR displacements from normal load (n2x, n3x, n3y):\n', np.hstack( (n2.v_disp[0], n3.v_disp) )

# Shear load
n3.F_ext = F_shear[1:]
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
print 'FEM displacements from shear load (n2x, n3x, n3y):\n', np.dot(Ki, F_shear)
print 'DR displacements from shear load (n2x, n3x, n3y):\n', np.hstack( (n2.v_disp[0], n3.v_disp) )

# Mixed load
n3.F_ext = F_mixed[1:]
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
print 'FEM displacements from mixed load (n2x, n3x, n3y):\n', np.dot(Ki, F_mixed)
print 'DR displacements from mixed load (n2x, n3x, n3y):\n', np.hstack( (n2.v_disp[0], n3.v_disp) ), '\n'

print '-------- nu = 0.2 --------\n\n'
e1.nu = 0.2
n3.F_ext = F_normal[1:]
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
K = e1.K[[2,4,5],:]
K = K[:,[2,4,5]]
Ki = np.linalg.inv(K)
print 'FEM displacements from normal load (n2x, n3x, n3y):\n', np.dot(Ki, F_normal)
print 'DR displacements from normal load (n2x, n3x, n3y):\n', np.hstack( (n2.v_disp[0], n3.v_disp) )

# Shear load
n3.F_ext = F_shear[1:]
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
print 'FEM displacements from shear load (n2x, n3x, n3y):\n', np.dot(Ki, F_shear)
print 'DR displacements from shear load (n2x, n3x, n3y):\n', np.hstack( (n2.v_disp[0], n3.v_disp) )

# Mixed load
n3.F_ext = F_mixed[1:]
dom.solve(dt=time_step, tol=tolerance, maxiter=maxiter, node_to_plot=node_to_plot)
print 'FEM displacements from mixed load (n2x, n3x, n3y):\n', np.dot(Ki, F_mixed)
print 'DR displacements from mixed load (n2x, n3x, n3y):\n', np.hstack( (n2.v_disp[0], n3.v_disp) ), '\n'
