'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., alfaC=0., stiffdef='isoparametric',
                             spring_stiffness_factor=10.,
                             load=[((49., 51.), (29., 31.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

dom.solve(dt=0.001, tol=1e-16, maxiter=2000, t_max=1., t_load=1., node_to_plot=59)
dom.plot(magnitude=200.)

print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[58].x, dom.nodes[58].y)
print 'CDEM DR deflection: %.6f' % dom.nodes[58].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'CDEM DR L2 norm: %.6f' % np.linalg.norm(v_disp)
print dom.get_stress_norm(ids=[0])
