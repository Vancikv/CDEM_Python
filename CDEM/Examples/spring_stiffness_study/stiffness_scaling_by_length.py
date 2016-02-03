'''
Created on Nov 24, 2015

@author: vancik
Same beam, different element sizes (>1.0), same contact stiffness.
Unscaled contacts: difference from FEM greater overall, greater for coarser mesh
Scaled contacts: difference from FEM lower overall, approximately the same for coarse and fine meshes

Alas, it does not really work that way
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np

for ni, nj in [[3, 10], [6, 20], [12, 40]]:
    Fz = -50.
    if nj != 10:
        Fz = -50. / (nj / 10 + 1)
    dom = get_2D_quadrangle_domain(ni=ni, nj=nj, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             domtype='CDEMStatic', spring_stiffness_factor=15., c_n=10000., c_s=10000.,
                             load=[((44.9, 55.1), (29., 31.), [0., Fz])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((99., 101.), (-1., 1.), [0, 1])])
    refnodes = dom.get_node_by_coords(50., 30., 1e-3)
    refnd = refnodes[0]
    dom.solve()
    # Unscaled
    u_defl = dom.nodes[refnd].v_disp[1]
    u_norm = dom.get_stress_norm([0])
    # Scaled
    dom.stiffness_scaling = True
    dom.solve()
    s_defl = dom.nodes[refnd].v_disp[1]
    s_norm = dom.get_stress_norm([0])
    # FEM
    domFEM = convert_CDEM_to_FEM_domain(dom)
    refnodes = domFEM.get_node_by_coords(50., 30., 1e-3)
    refnd = refnodes[0]
    for nd in domFEM.nodes:
        nd.F_ext[0] *= 2.
        nd.F_ext[1] *= 2.
    domFEM.solve()
    F_defl = domFEM.nodes[refnd].v_disp[1]
    F_norm = domFEM.get_stress_norm([0])

    print 'Mesh %d x %d:' % (ni, nj)
    print '                    Unscaled (/FEM)             Scaled (/FEM)             FEM'
    print '    Deflections:    %.6f (%.6f)        %.6f (%.6f)      %.6f' % (u_defl, u_defl / F_defl, s_defl, s_defl / F_defl, F_defl)
    print '    Norms:          %.6f (%.6f)        %.6f (%.6f)      %.6f' % (u_norm, u_norm / F_norm, s_norm, s_norm / F_norm, F_norm)


