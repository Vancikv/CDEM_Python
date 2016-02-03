'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *

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
