'''
Created on 11. 11. 2015

@author: Kapsak
'''

    def set_K_isoparametric(self):
        '''
        Reduced integration with hourglass stabilization is used, B_red = j_0_dev * L_hg * M_hg
        '''
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)

        # Get node info
        x = np.array([self.domain().nodes[i - 1].x for i in self.nodes])
        y = np.array([self.domain().nodes[i - 1].y for i in self.nodes])

        # Get Jacobi matrix at center point
        g_xi = np.array([-0.25, 0.25, 0.25, -0.25])
        g_eta = np.array([-0.25, -0.25, 0.25, 0.25])
        h = np.array([.25, -.25, .25, -.25])
        J_0 = np.array([[np.dot(g_xi, x), np.dot(g_eta, x)],
                        [np.dot(g_xi, y), np.dot(g_eta, y)]])
        J_0_inv = np.linalg.inv(J_0)
        det_J_0 = np.linalg.det(J_0)

        # Get Jacobi matrices at gauss points
        gp = (1. / 3.) ** 0.5
        gps = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
        J = [np.array([[np.dot(g_xi + eta * h, x), np.dot(g_eta + xi * h, x)],
                        [np.dot(g_xi + eta * h, y), np.dot(g_eta + xi * h, y)]]) for xi, eta in gps]
        J_inv = [np.linalg.inv(m) for m in J]

        # Get B_0 matrix
        b = np.dot(np.transpose(np.vstack((g_xi, g_eta))), J_0_inv)
        B_0 = np.vstack((np.hstack((b[:, 0], np.zeros(4))), np.hstack((np.zeros(4), b[:, 1])), np.hstack((b[:, 1], b[:, 0]))))

        # Get B_hg matrix
        gamma = (np.eye(4) - np.dot(np.dot(np.transpose(np.vstack((g_xi, g_eta))), J_0_inv), np.vstack((x, y))))
        gamma = np.dot(gamma, h)
        j_0_dev = 1. / 3. * np.vstack((np.hstack((2 * J_0_inv[:, 0], -1 * J_0_inv[:, 1])),
                             np.hstack((-1 * J_0_inv[:, 0], 2 * J_0_inv[:, 1])),
                             np.hstack((3 * J_0_inv[:, 0], 3 * J_0_inv[:, 1]))))
        L_hg = [np.array([[eta, 0.], [xi, 0.], [0., eta], [0., xi]]) for xi, eta in gps]
        M_hg = np.vstack((np.hstack((gamma, np.zeros(4))), np.hstack((np.zeros(4), gamma))))
        B_red = [np.dot(np.dot(j_0_dev, mL), M_hg) for mL in L_hg]
        self.B = [B_0 + b_red for b_red in B_red]
        for i in range(len(self.B)):
            self.B[i] = self.B[i][:, [0, 4, 1, 5, 2, 6, 3, 7]]

        # Stiffness matrix
        self.volume = 4 * det_J_0 * self.thickness
        K_red = self.volume / 4. * sum(np.dot(np.dot(np.transpose(B), C), B) for B in B_red)
        K = K_red + self.volume * np.dot(np.dot(np.transpose(B_0), C) , B_0)

        # Rearrange to fit vector [ux1,uy1,ux2,uy2,ux3,uy3,ux4,uy4]
        K = K[[0, 4, 1, 5, 2, 6, 3, 7], :]
        K = K[:, [0, 4, 1, 5, 2, 6, 3, 7]]
        self.K = K

    def set_K_isoparametric(self):
        '''
        Full integration of hourglass part => shear locking
        '''
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)

        # Get node info
        x = np.array([self.domain().nodes[i - 1].x for i in self.nodes])
        y = np.array([self.domain().nodes[i - 1].y for i in self.nodes])

        # Get Jacobi matrix at center point
        g_xi = np.array([-0.25, 0.25, 0.25, -0.25])
        g_eta = np.array([-0.25, -0.25, 0.25, 0.25])
        h = np.array([.25, -.25, .25, -.25])
        J_0 = np.array([[np.dot(g_xi, x), np.dot(g_eta, x)],
                        [np.dot(g_xi, y), np.dot(g_eta, y)]])
        J_0_inv = np.linalg.inv(J_0)
        det_J_0 = np.linalg.det(J_0)

        # Get Jacobi matrices at gauss points
        gp = (1. / 3.) ** 0.5
        gps = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
        J = [np.array([[np.dot(g_xi + eta * h, x), np.dot(g_eta + xi * h, x)],
                        [np.dot(g_xi + eta * h, y), np.dot(g_eta + xi * h, y)]]) for xi, eta in gps]
        J_inv = [np.linalg.inv(m) for m in J]

        # Get B_0 matrix
        b = np.dot(np.transpose(np.vstack((g_xi, g_eta))), J_0_inv)
        B_0 = np.vstack((np.hstack((b[:, 0], np.zeros(4))), np.hstack((np.zeros(4), b[:, 1])), np.hstack((b[:, 1], b[:, 0]))))

        # Get B_hg matrix
        gamma = (np.eye(4) - np.dot(np.dot(np.transpose(np.vstack((g_xi, g_eta))), J_0_inv), np.vstack((x, y))))
        gamma = np.dot(gamma, h)
        j = [np.vstack((np.hstack((m[0, :], np.zeros(2))), np.hstack((np.zeros(2), m[1, :])), np.hstack((m[1, :], m[0, :])))) for m in J_inv]
        L_hg = [np.array([[eta, 0.], [xi, 0.], [0., eta], [0., xi]]) for xi, eta in gps]
        M_hg = np.vstack((np.hstack((gamma, np.zeros(4))), np.hstack((np.zeros(4), gamma))))
        B_hg = [np.dot(np.dot(mj, mL), M_hg) for mj, mL in zip(j, L_hg)]

        # Stiffness matrix
        self.volume = 4 * det_J_0 * self.thickness
        K_hg = self.volume / 4. * sum(np.dot(np.dot(np.transpose(B), C), B) for B in B_hg)
        K = self.volume * np.dot(np.dot(np.transpose(B_0), C) , B_0) + K_hg

        # Rearrange to fit vector [ux1,uy1,ux2,uy2,ux3,uy3,ux4,uy4]
        K = K[[0, 4, 1, 5, 2, 6, 3, 7], :]
        K = K[:, [0, 4, 1, 5, 2, 6, 3, 7]]
        self.K = K

    def set_K_isoparametric(self):
        '''
        Reduced integration with hourglass stabilization is used, B_red = j_0_dev * L_hg * M_hg
        '''
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)

        # Get node info
        x = np.array([self.domain().nodes[i - 1].x for i in self.nodes])
        y = np.array([self.domain().nodes[i - 1].y for i in self.nodes])

        # Get Jacobi matrix at center point
        g_xi = np.array([-0.25, 0.25, 0.25, -0.25])
        g_eta = np.array([-0.25, -0.25, 0.25, 0.25])
        h = np.array([.25, -.25, .25, -.25])
        J_0 = np.array([[np.dot(g_xi, x), np.dot(g_eta, x)],
                        [np.dot(g_xi, y), np.dot(g_eta, y)]])
        J_0_inv = np.linalg.inv(J_0)
        det_J_0 = np.linalg.det(J_0)

        # Get Jacobi matrices at gauss points
        gp = (1. / 3.) ** 0.5
        gps = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
        J = [np.array([[np.dot(g_xi + eta * h, x), np.dot(g_eta + xi * h, x)],
                        [np.dot(g_xi + eta * h, y), np.dot(g_eta + xi * h, y)]]) for xi, eta in gps]
        J_inv = [np.linalg.inv(m) for m in J]

        # Get B_0 matrix
        b = np.dot(np.transpose(np.vstack((g_xi, g_eta))), J_0_inv)
        B_0 = np.vstack((np.hstack((b[:, 0], np.zeros(4))), np.hstack((np.zeros(4), b[:, 1])), np.hstack((b[:, 1], b[:, 0]))))

        # Get B_hg matrix
        gamma = (np.eye(4) - np.dot(np.dot(np.transpose(np.vstack((g_xi, g_eta))), J_0_inv), np.vstack((x, y))))
        gamma = np.dot(gamma, h)
        j_0_dev = 1. / 3. * np.vstack((np.hstack((2 * J_0_inv[:, 0], -1 * J_0_inv[:, 1])),
                             np.hstack((-1 * J_0_inv[:, 0], 2 * J_0_inv[:, 1])),
                             np.hstack((3 * J_0_inv[:, 0], 3 * J_0_inv[:, 1]))))
        L_hg = [np.array([[eta, 0.], [xi, 0.], [0., eta], [0., xi]]) for xi, eta in gps]
        M_hg = np.vstack((np.hstack((gamma, np.zeros(4))), np.hstack((np.zeros(4), gamma))))
        B_red = [np.dot(np.dot(j_0_dev, mL), M_hg) for mL in L_hg]
        self.B = [B_0 + b_red for b_red in B_red]

        # Stiffness matrix
        self.volume = 4 * det_J_0 * self.thickness
        K_red = self.volume / 4. * sum(np.dot(np.dot(np.transpose(B), C), B) for B in B_red)
        K = self.volume * np.dot(np.dot(np.transpose(B_0), C) , B_0) + K_red

        # Rearrange to fit vector [ux1,uy1,ux2,uy2,ux3,uy3,ux4,uy4]
        K = K[[0, 4, 1, 5, 2, 6, 3, 7], :]
        K = K[:, [0, 4, 1, 5, 2, 6, 3, 7]]
        self.K = K

'''
Created on Nov 24, 2015

@author: vancik
Same beam, different element sizes (>1.0), same contact stiffness.
Unscaled contacts: difference from FEM greater overall, greater for coarser mesh
Scaled contacts: difference from FEM lower overall, approximately the same for coarse and fine meshes
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np

domf = get_2D_quadrangle_domain(ni=12, nj=40, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=3., domtype='CDEMStatic',
                             load=[((49., 51.), (29., 31.), [0., -50.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((99., 101.), (-1., 1.), [0, 1])])
domc = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             spring_stiffness_factor=3., domtype='CDEMStatic',
                             load=[((49., 51.), (29., 31.), [0., -50.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((99., 101.), (-1., 1.), [0, 1])])
domfFEM = convert_CDEM_to_FEM_domain(domf)
refnodes = domfFEM.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
domfFEM.nodes[refnd].F_ext = [0.,-100.]
domcFEM = convert_CDEM_to_FEM_domain(domc)
refnodes = domcFEM.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
domcFEM.nodes[refnd].F_ext = [0.,-100.]

for ni,nj in [[3,10],[6,20],[12,40]]:
    
print 'UNSCALED FINE:\n'
domf.solve()
refnodes = domf.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
#print 'Reference node coordinates: x = %f, y = %f' % (domf.nodes[refnd].x, domf.nodes[refnd].y)
fu_defl = domf.nodes[refnd].v_disp[1]
fu_norm = domf.get_stress_norm([0])
print 'CDEM static fine mesh deflection: %.6f' % fu_defl
print 'CDEM static fine mesh sigma_x L2 norm: %.6f' % fu_norm
#domf.plot(magnitude=50.)

print 'SCALED FINE:\n'
domf.stiffness_scaling = True
domf.solve()
#print 'Reference node coordinates: x = %f, y = %f' % (domf.nodes[refnd].x, domf.nodes[refnd].y)
fs_defl = domf.nodes[refnd].v_disp[1]
fs_norm = domf.get_stress_norm([0])
print 'CDEM static fine mesh deflection: %.6f' % fs_defl
print 'CDEM static fine mesh sigma_x L2 norm: %.6f' % fs_norm
#domf.plot(magnitude=50.)

print 'FEM FINE:\n'
domfFEM.solve()
refnodes = domfFEM.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
#print 'Reference node coordinates: x = %f, y = %f' % (domfFEM.nodes[refnd].x, domfFEM.nodes[refnd].y)
fF_defl = domfFEM.nodes[refnd].v_disp[1]
fF_norm = domfFEM.get_stress_norm([0])
print 'FEM fine mesh deflection: %.6f' % fF_defl
print 'FEM fine mesh sigma_x L2 norm: %.6f' % fF_norm
#domfFEM.plot(magnitude=50.)
print '\nUnscaled deflection ratio: %.6f' % (fu_defl / fF_defl)
print 'Unscaled norm ratio: %.6f' % (fu_norm / fF_norm)
print 'Scaled deflection ratio: %.6f' % (fs_defl / fF_defl)
print 'Scaled norm ratio: %.6f' % (fs_norm / fF_norm)


print '\n\n\nUNSCALED COARSE:\n'
domc.solve()
refnodes = domc.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
#print 'Reference node coordinates: x = %f, y = %f' % (domc.nodes[refnd].x, domc.nodes[refnd].y)
cu_defl = domc.nodes[refnd].v_disp[1]
cu_norm = domc.get_stress_norm([0])
print 'CDEM static coarse mesh deflection: %.6f' % cu_defl
print 'CDEM static coarse mesh sigma_x L2 norm: %.6f' % cu_norm
#domc.plot(magnitude=50.)

print 'SCALED COARSE:\n'
domc.stiffness_scaling = True
domc.solve()
#print 'Reference node coordinates: x = %f, y = %f' % (domc.nodes[refnd].x, domc.nodes[refnd].y)
cs_defl = domc.nodes[refnd].v_disp[1]
cs_norm = domc.get_stress_norm([0])
print 'CDEM static coarse mesh deflection: %.6f' % cs_defl
print 'CDEM static coarse mesh sigma_x L2 norm: %.6f' % cs_norm
#domc.plot(magnitude=50.)

print 'FEM COARSE:\n'
domcFEM.solve()
refnodes = domcFEM.get_node_by_coords(50., 30., 1e-3)
refnd = refnodes[0]
#print 'Reference node coordinates: x = %f, y = %f' % (domcFEM.nodes[refnd].x, domcFEM.nodes[refnd].y)
cF_defl = domcFEM.nodes[refnd].v_disp[1]
cF_norm = domcFEM.get_stress_norm([0])
print 'FEM coarse mesh deflection: %.6f' % cF_defl
print 'FEM coarse mesh sigma_x L2 norm: %.6f' % cF_norm
#domcFEM.plot(magnitude=50.)

print '\nUnscaled deflection ratio: %.6f' % (cu_defl / cF_defl)
print 'Unscaled norm ratio: %.6f' % (cu_norm / cF_norm)
print 'Scaled deflection ratio: %.6f' % (cs_defl / cF_defl)
print 'Scaled norm ratio: %.6f' % (cs_norm / cF_norm)

