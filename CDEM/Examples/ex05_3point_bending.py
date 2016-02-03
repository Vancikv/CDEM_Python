'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.0,
                             thickness=1., alfaC=100., stiffdef='isoparametric',
                             c_n=50000., c_s=50000.,
                             load=[((45., 55.), (29., 31.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

dom.solve(dt=0.01, tol=1e-16, maxiter=1000, node_to_plot=59)
dom.plot(magnitude=100.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[58].x, dom.nodes[58].y)
print 'CDEM DR deflection: %.6f' % dom.nodes[58].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'CDEM DR L2 norm: %.6f' % np.linalg.norm(v_disp)
print 'CDEM DR sigma_y norm: %.6f' % dom.get_stress_norm(ids=[0])

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             c_n=50000., c_s=50000., domtype='CDEMStatic',
                             load=[((45., 55.), (29., 31.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

dom.solve()
dom.plot(magnitude=100.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[58].x, dom.nodes[58].y)
print 'CDEM static deflection: %.6f' % dom.nodes[58].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'CDEM static L2 norm: %.6f' % np.linalg.norm(v_disp)
print 'CDEM static sigma_y norm: %.6f' % dom.get_stress_norm(ids=[0])

dom = get_2D_quadrangle_FEM_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             load=[((45., 55.), (29., 31.), [0., -100.])],
                             supports=[((-5., 5.), (-5., 5.), [1, 1]), ((95., 105.), (-5., 5.), [0, 1])])

dom.solve()
dom.plot(magnitude=100.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[23].x, dom.nodes[23].y)
print 'FEM static deflection: %.6f' % dom.nodes[23].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'FEM static L2 norm: %.6f' % np.linalg.norm(v_disp)
print 'FEM static sigma_y norm: %.6f' % dom.get_stress_norm(ids=[0])
'''
15k steps:
Maximum step count reached: maximum criterium = 0.0000027516
Reference node coordinates: x = 50.000000, y = 30.000000
CDEM DR deflection: -0.205491
CDEM DR L2 norm: 1.577401
CDEM DR sigma_y norm: 89.0520080706
Reference node coordinates: x = 50.000000, y = 30.000000
CDEM static deflection: -0.145975
CDEM static L2 norm: 1.203950
CDEM static sigma_y norm: 59.467031
Reference node coordinates: x = 50.000000, y = 30.000000
FEM static deflection: -0.132265
FEM static L2 norm: 0.646676
FEM static sigma_y norm: 59.732348

50k steps:
Maximum step count reached: maximum criterium = 0.0000002034
Reference node coordinates: x = 50.000000, y = 30.000000
CDEM DR deflection: -0.342450
CDEM DR L2 norm: 2.945757
CDEM DR sigma_y norm: 147.259735

100k steps:
Maximum step count reached: maximum criterium = 0.0000000061
Reference node coordinates: x = 50.000000, y = 30.000000
CDEM DR deflection: -0.385888
CDEM DR L2 norm: 3.414569
CDEM DR sigma_y norm: 167.454235
'''