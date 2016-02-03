'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM import *

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=10., domtype='CDEMStatic',
                             load=[((49, 51.), (29., 31.), [0., -50.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((99., 101.), (-1., 1.), [0, 1])])

dom.solve()
dom.plot(magnitude=5.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[58].x, dom.nodes[58].y)
print 'CDEM static deflection: %.6f' % dom.nodes[58].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'CDEM static L2 norm: %.6f' % np.linalg.norm(v_disp)
print dom.get_stress_norm(ids=[0])
