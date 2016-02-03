'''
Created on Jan 4, 2016

@author: vancik
'''

import numpy as np
import matplotlib.pyplot as plt
import time
from basic import *

class Element(object):
    def __init__(self, domain=None, nodes=[], E=0., nu=0., density=0., thickness=1., alfaC=1., stiffdef=None, *args, **kwargs):
        # Nodes numbered from 1
        self.nodes = nodes
        self.E = E
        self.nu = nu
        self.density = density
        self.thickness = thickness
        self.alfaC = alfaC
        self.domain = domain
        if stiffdef:
            self.stiffdef = stiffdef

    def set_matrices(self):
        '''
        Calculate and store local stiffness, mass and damping matrices
        '''
        pass

    def save_next_vals(self, v_disp, v_velo, v_acce):
        self.v_disp = v_disp
        self.v_velo = v_velo
        self.v_acce = v_acce

    def write_next_vals(self):
        dom = self.domain()
        nds = [dom.nodes[i - 1] for i in self.nodes]
        for i, nd in enumerate(nds):
            nd.set_disp(self.v_disp[2 * i:2 * i + 2])
            nd.set_velo(self.v_velo[2 * i:2 * i + 2])
            nd.set_acce(self.v_acce[2 * i:2 * i + 2])
        pass

    def calc_normal_vectors(self):
        '''
        Calculate the normal vector of each face and pass it to the adjacent nodes
        '''
        dom = self.domain()
        nds = self.nodes + [self.nodes[0]]
        x = [dom.nodes[i - 1].x for i in nds]
        y = [dom.nodes[i - 1].y for i in nds]
        for i in range(len(nds) - 1):
            dx, dy = x[i + 1] - x[i], y[i + 1] - y[i]
            # Unit normal vector pointing outwards when the nodes are counter-clockwise
            L = (dx * dx + dy * dy) ** 0.5
            n = 1 / L * np.array([dy, -dx])
            dom.nodes[nds[i] - 1].n_vects[1] = n
            dom.nodes[nds[i] - 1].l_edges[1] = L
            dom.nodes[nds[i + 1] - 1].n_vects[0] = n
            dom.nodes[nds[i + 1] - 1].l_edges[0] = L

    def set_codes(self):
        self.v_code = reduce(lambda x, y: x + y, [self.domain().nodes[i - 1].v_code for i in self.nodes])

    def plot(self, ax, magnitude):
        dom = self.domain()
        nds = self.nodes + [self.nodes[0]]
        xy = np.transpose(np.array([[dom.nodes[i - 1].x, dom.nodes[i - 1].y] for i in nds]))
        u = np.transpose(np.array([dom.nodes[i - 1].v_disp for i in nds]))

        # Plot reference shape
        ax.plot(xy[0, :], xy[1, :], linestyle='-', color='black', linewidth=1.)
        # Plot deformed shape
        ax.plot(xy[0, :] + magnitude * u[0, :], xy[1, :] + magnitude * u[1, :], linestyle='-', color='blue', linewidth=3.)

    def iterate(self, dt, tau, timelog=None, verbose=False):
        dom = self.domain()
        nds = [dom.nodes[i - 1] for i in self.nodes]

        t1 = time.clock()
        # Obtain nodal values
        v_disp = np.hstack([nd.v_disp for nd in nds])  # Displacement vector
        v_velo = np.hstack([nd.v_velo for nd in nds])  # Velocity vector
        v_acce = np.hstack([nd.v_acce for nd in nds])  # Acceleration vector

        t2 = time.clock()

        # Calc next values
        v_disp = v_disp + dt * v_velo + 0.5 * dt * dt * v_acce
        v_velo = v_velo + dt * v_acce

        t3 = time.clock()
        # Calc forces
        F_k_e = -1 * np.dot(self.K, v_disp)  # Element stiffness forces
        F_r = -1 * F_k_e * np.hstack([nd.supports for nd in nds])  # Reaction forces
        F_k_c = np.hstack([dom.get_contact_force(i) for i in self.nodes])  # Contact forces
        F_c = -1. * np.dot(self.C, v_velo)  # Damping forces
        F_g = load_function(tau) * np.hstack([nd.F_g for nd in nds])  # Gravity forces
        F_ext = load_function(tau) * np.hstack([nd.F_ext for nd in nds])  # External forces
        F_tot = F_k_e + F_k_c + F_c + F_g + F_ext + F_r  # Resulting forces

        if verbose:
            print 'Displacement', v_disp
            print 'Velocity', v_velo
            print 'Element stiffness\n', F_k_e
            print 'Contact stiffness\n', F_k_c
            print 'Damping\n', F_c, '\n'

        t4 = time.clock()
        # Calc next acceleration
        v_acce = np.dot(self.M_inv, F_tot)

        t5 = time.clock()
        # Save next values
        self.save_next_vals(v_disp, v_velo, v_acce)
        t6 = time.clock()
        if timelog:
            timelog.write('    Obtaining nodal values: %g\n' % (t2 - t1))
            timelog.write('    Calculating next values: %g\n' % (t3 - t2))
            timelog.write('    Calculating forces: %g\n' % (t4 - t3))
            timelog.write('    Calculating next acceleration: %g\n' % (t5 - t4))
            timelog.write('    Writing next values: %g\n' % (t6 - t5))

        return np.dot(v_velo, v_velo)  # Convergence criterium

    def iterate_kin_damping(self, dt, tau, verbose=False):
        dom = self.domain()
        v_disp = np.hstack([dom.nodes[i - 1].v_disp for i in self.nodes])  # Displacement vector
        v_velo = np.hstack([dom.nodes[i - 1].v_velo for i in self.nodes])  # Velocity vector
        v_acce = np.hstack([dom.nodes[i - 1].v_acce for i in self.nodes])  # Acceleration vector
        if verbose:
            print 'Displacement', v_disp
            print 'Velocity', v_velo
            print 'Acceleration', v_acce
        v_disp = v_disp + dt * v_velo + 0.5 * dt * dt * v_acce
        v_velo = v_velo + dt * v_acce
        F_k_e = -1 * np.dot(self.K, v_disp)  # Element stiffness forces
        F_r = -1 * F_k_e * np.hstack([dom.nodes[i - 1].supports for i in self.nodes])  # Reaction forces
        F_k_c = np.hstack([dom.get_contact_force(i) for i in self.nodes])  # Contact forces
        if verbose:
            print 'Displacement', v_disp
            print 'Velocity', v_velo
            print 'Acceleration', v_acce
            # print 'Element stiffness\n', F_k_e
            # print 'Contact stiffness\n', F_k_c
        F_g = load_function(tau) * np.hstack([dom.nodes[i - 1].F_g for i in self.nodes])  # Gravity forces
        F_ext = load_function(tau) * np.hstack([dom.nodes[i - 1].F_ext for i in self.nodes])  # External forces

        F_tot = F_k_e + F_k_c + F_g + F_ext + F_r  # Resulting forces
        v_acce = np.dot(np.linalg.inv(self.M), F_tot)

        return (np.dot(v_velo, v_velo), v_disp, v_velo, v_acce)  # kinetic energy

    def set_node_vals(self, disp, velo, acce):
        dom = self.domain()
        for i, j in enumerate(self.nodes):
            dom.nodes[j - 1].set_disp(disp[2 * i:2 * i + 2])
            dom.nodes[j - 1].set_velo(velo[2 * i:2 * i + 2])
            dom.nodes[j - 1].set_acce(acce[2 * i:2 * i + 2])

class ElemTriangleLin(Element):
    def set_matrices(self, verbose=False):
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)

        # Get node info
        x = np.array([self.domain().nodes[i - 1].x for i in self.nodes])
        y = np.array([self.domain().nodes[i - 1].y for i in self.nodes])

        rot = lambda l, n: l[n:] + l[:n]
        ids = [0, 1, 2]
        a = [x[j] * y[k] - x[k] * y[j] for i, j, k in [rot(ids, n) for n in ids]]
        b = [y[j] - y[k] for i, j, k in [rot(ids, n) for n in ids]]
        c = [x[k] - x[j] for i, j, k in [rot(ids, n) for n in ids]]
        self.A = 0.5 * np.linalg.det(np.transpose(np.vstack((np.ones(3), x, y))))
        A = self.A

        # du/dx, du/dy, B matrix
        b1 = 1 / 2. / A * np.array(b)
        b2 = 1 / 2. / A * np.array(c)
        B = np.vstack((np.hstack((b1, np.zeros(3))), np.hstack((np.zeros(3), b2)), np.hstack((b2, b1))))
        self.B = [B]

        K = A * self.thickness * np.dot(np.dot(np.transpose(B), C) , B)
        K = K[[0, 3, 1, 4, 2, 5], :]
        K = K[:, [0, 3, 1, 4, 2, 5]]
        self.K = K
        self.M = (A * self.thickness * self.density / 3.) * np.eye(6, 6)
        self.M_inv = 1. / (A * self.thickness * self.density / 3.) * np.eye(6, 6)
        for i in range(len(self.nodes)):
            self.domain().nodes[self.nodes[i] - 1].mass = self.M[i, i]
        # Damping matrix
        # self.C = self.M * self.alfaC
        self.C = self.alfaC * np.diag([(K[i, i] / self.M[i, i]) ** 0.5 for i in range(K.shape[0])])
        if verbose: print self.K, self.M, self.C

    def get_stress(self):
        dom = self.domain()
        v_disp = np.hstack([dom.nodes[i - 1].v_disp for i in self.nodes])
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)
        v_sigma_arr = [np.dot(C, np.dot(Bi, v_disp)) for Bi in self.B]
        return v_sigma_arr

class ElemQuadrangleLin(Element):
    stiffdef = 'isoparametric'
    stiffness_matrix_definition = {'2triangles':'set_K_2_triangles', '4triangles':'set_K_4_triangles',
                                   'isoparametric':'set_K_isoparametric'}

    def set_K_2_triangles(self):
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)

        nds = self.nodes + [self.nodes[0]]
        x = np.array([self.domain().nodes[i - 1].x for i in nds])
        y = np.array([self.domain().nodes[i - 1].y for i in nds])

        data = []
        data.append(get_triangle_stiffness_matrix(C, x[0:3], y[0:3], self.thickness))
        data.append(get_triangle_stiffness_matrix(C, x[[2, 3, 0]], y[[2, 3, 0]], self.thickness))
        k_lst, A_lst = zip(*data)

        K = np.zeros((8, 8))
        K[:6, :6] += k_lst[0]

        K[4:, 4:] += k_lst[1][:4, :4]
        K[:2, :2] += k_lst[1][4:, 4:]
        K[:2, 4:] += k_lst[1][4:, :4]
        K[4:, :2] += k_lst[1][:4, 4:]
        self.K = K
        self.volume = self.thickness * sum(A_lst)

    def set_K_4_triangles(self):
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)

        nds = self.nodes
        x = np.array([self.domain().nodes[i - 1].x for i in nds])
        y = np.array([self.domain().nodes[i - 1].y for i in nds])

        # Auxiliary center node
        xs = .25 * sum(x)
        ys = .25 * sum(y)

        # Triangle nodes
        xt_lst = [np.hstack((x[[i, (i + 1) % 4]], xs)) for i in range(4)]
        yt_lst = [np.hstack((y[[i, (i + 1) % 4]], xs)) for i in range(4)]

        data = [get_triangle_stiffness_matrix(C, xt, yt, self.thickness) for xt, yt in zip(xt_lst, yt_lst)]
        k_lst, A_lst = zip(*data)

        K = np.zeros((10, 10))
        for i in range(4):
            ids = [2 * i, 2 * i + 1, (2 * i + 2) % 8, (2 * i + 3) % 8, 8, 9]
            for jj, j in enumerate(ids):
                for kk, k in enumerate(ids):
                    K[j, k] += k_lst[i][jj, kk]
        for i in range(8):
            K[i, :] -= K[i, 9] / K[9, 9] * K[9, :] + K[i, 8] / K[8, 8] * K[8, :]
        self.K = K[:8, :8]
        self.volume = sum(A_lst) * self.thickness

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

    def get_stress(self):
        dom = self.domain()
        v_disp = np.hstack([dom.nodes[i - 1].v_disp for i in self.nodes])
        C = get_plane_stress_stiffness_matrix(self.E, self.nu)
        v_sigma_arr = [np.dot(C, np.dot(Bi, v_disp)) for Bi in self.B]
        return v_sigma_arr

    def set_matrices(self, verbose=False):
        '''
        Calculate and store local stiffness, mass and damping matrices.
        '''
        getattr(self, self.stiffness_matrix_definition[self.stiffdef])()
        # Mass matrix + set node masses
        self.M = (self.volume * self.density / 4.) * np.eye(8, 8)
        self.M_inv = 1. / (self.volume * self.density / 4.) * np.eye(8, 8)
        for i in range(len(self.nodes)):
            self.domain().nodes[self.nodes[i] - 1].mass = self.M[i, i]
        # Damping matrix
        # self.C = self.M * self.alfaC
        self.C = self.alfaC * np.diag([(self.K[i, i] / self.M[i, i]) ** 0.5 for i in range(self.K.shape[0])])
        if verbose: print 'Element stiffness matrix:\n', self.K, 'Element mass matrix:\n', self.M, 'Element damping matrix:\n', self.C
