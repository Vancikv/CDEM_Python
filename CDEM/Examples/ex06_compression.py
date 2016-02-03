'''
Created on 22. 11. 2015

@author: Kapsak
'''


from CDEM.basic import *
from CDEM.utils import *

dom = get_2D_quadrangle_domain(ni=10, nj=4, lx=20., ly=100., E=10000., nu=0.2, density=1.5,
                             thickness=1., alfaC=100., stiffdef='isoparametric',
                             c_n=50000., c_s=50000.,
                             load=[((-1., 21.), (99., 101.), [0., -125.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((1., 21.), (-1., 1.), [0, 1])])

dom.solve(dt=0.01, tol=1e-16, maxiter=5, node_to_plot=78)
dom.plot(magnitude=10.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[78].x, dom.nodes[78].y)
print 'CDEM DR deflection: %.6f' % dom.nodes[78].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'CDEM DR L2 norm: %.6f' % np.linalg.norm(v_disp)
print 'CDEM DR sigma_y norm: %.6f' % dom.get_stress_norm(ids=[1])

dom = get_2D_quadrangle_domain(ni=10, nj=4, lx=20., ly=100., E=10000., nu=0.2, density=1.,
                             thickness=1., alfaC=100., stiffdef='isoparametric',
                             c_n=50000., c_s=50000., domtype='CDEMStatic',
                             load=[((-1., 21.), (99., 101.), [0., -125.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((1., 21.), (-1., 1.), [0, 1])])

dom.solve()
dom.plot(magnitude=10.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[78].x, dom.nodes[78].y)
print 'CDEM static deflection: %.6f' % dom.nodes[78].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'CDEM static L2 norm: %.6f' % np.linalg.norm(v_disp)
print 'CDEM static sigma_y norm: %.6f' % dom.get_stress_norm(ids=[1])

dom = get_2D_quadrangle_FEM_domain(ni=10, nj=4, lx=20., ly=100., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             load=[((-1., 21.), (99., 101.), [0., -200.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((1., 21.), (-1., 1.), [0, 1])])

dom.solve()
dom.plot(magnitude=10.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[21].x, dom.nodes[21].y)
print 'FEM static deflection: %.6f' % dom.nodes[21].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'FEM static L2 norm: %.6f' % np.linalg.norm(v_disp)
print 'FEM static sigma_y norm: %.6f' % dom.get_stress_norm(ids=[1])
'''
15k steps:
Maximum step count reached: maximum criterium = 0.0000010028
Reference node coordinates: x = 10.000000, y = 100.000000
CDEM DR deflection: -0.508889
CDEM DR L2 norm: 3.690920
CDEM DR sigma_y norm: 601.057014
Reference node coordinates: x = 10.000000, y = 100.000000
CDEM static deflection: -0.536054
CDEM static L2 norm: 3.949661
CDEM static sigma_y norm: 632.740232

50k steps:
Maximum step count reached: maximum criterium = 0.0000000007
Reference node coordinates: x = 10.000000, y = 100.000000
CDEM DR deflection: -0.535950
CDEM DR L2 norm: 3.930439
CDEM DR sigma_y norm: 632.732874
Reference node coordinates: x = 10.000000, y = 100.000000
CDEM static deflection: -0.536054
CDEM static L2 norm: 3.949661
CDEM static sigma_y norm: 632.740232

'''

