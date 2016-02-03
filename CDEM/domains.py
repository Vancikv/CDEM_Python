'''
Created on Jan 4, 2016

@author: vancik
'''

import numpy as np
import weakref
import matplotlib.pyplot as plt
import itertools as it
import tabulate as tbl
import time
from basic import load_function

class Domain(object):
    def __init__(self, elements=[], nodes=[]):
        self.elements = elements
        for el in self.elements:
            el.domain = weakref.ref(self)
        self.nodes = nodes

    def plot(self, magnitude=10., savefile='', show=True):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for el in self.elements:
            el.plot(ax, magnitude)
        plt.axis('equal')
        if savefile:
            fig.savefig(savefile)
        if show:
            plt.show()

    def save_to_text_file(self, path):
        file = open(path, 'w')
        nds = self.nodes
        els = self.elements
        file.write("nnodes %d nelems %d cn %g cs %g\n" % (len(nds), len(self.elements), self.C[0, 0], self.C[1, 1]))
        for i, nd in enumerate(nds):
            file.write("node %d ndofs %d position %g %g neighbors %d %d supports %d %d load %g %g\n"
                       % (i + 1, nd.nnodedofs, nd.x, nd.y, nd.neighbors[0], nd.neighbors[1], nd.supports[0], nd.supports[1], nd.F_ext[0], nd.F_ext[1]))
        for i, el in enumerate(els):
            file.write("element %d nodes " % (i + 1))
            file.write("%d " % len(el.nodes))
            for j in range(len(el.nodes)):
                file.write("%d " % (el.nodes[j]))
            file.write("E %g nu %g density %g thickness %g alfaC %g\n" % (el.E, el.nu, el.density, el.thickness, el.alfaC))

    def get_node_by_coords(self, x, y, rad=1e-6):
        '''
        Return a list of node numbers within the specified radius from point x,y
        '''
        nds = []
        isclose = lambda n1, n2, tol = 1e-6: (((n2[0] - n1[0]) ** 2 + (n2[1] - n1[1]) ** 2) ** 0.5) < tol
        for i, nd in enumerate(self.nodes):
            if isclose((x, y), (nd.x, nd.y), tol=rad):
                nds.append(i)
        return nds

    def solve(self):
        pass

    def get_stress_norm(self, ids=[0, 1, 2]):
        '''
        Calculates the L2 norm of chosen stress components across all the elements.
        ids is a list of stress component ids: 0 - sigma_x, 1 - sigma_y, 2 - tau_xy
        '''
        norm = 0.
        for el in self.elements:
            v_sigma_arr = el.get_stress()
            for v_sigma in v_sigma_arr:
                norm += np.dot(v_sigma[ids], v_sigma[ids])
        return norm ** 0.5

class DomainCDEMDR(Domain):
    def __init__(self, c_n, c_s, contact_stiffness='uniform', elements=[], nodes=[]):
        self.elements = elements
        for el in self.elements:
            el.domain = weakref.ref(self)
        self.nodes = nodes
        self.set_contacts(c_n, c_s)

    def set_contacts(self, c_n, c_s):
        self.C = np.array([[c_n, 0.],
                           [0., c_s]])  # Contact stiffness matrix

    def get_contact_force(self, node):
        n0 = self.nodes[node - 1]
        F = np.array([0., 0.])
        for nbr, n_vect in zip(n0.neighbors, n0.n_vects):
            if nbr:
                T = np.array([[n_vect[0], n_vect[1]],
                              [-n_vect[1], n_vect[0]]])  # Transformation matrix
                du_g = self.nodes[nbr - 1].v_disp - n0.v_disp
                F += np.dot(np.dot(np.transpose(T), self.C), np.dot(T, du_g))
        return F

    def solve(self, dt=0.01, tol=0.001, maxiter=50, t_max=2., t_load=1., node_to_plot=None, verbose=False,
              savefile='', savetextfile='', show=True, timelog=None, output_freq=1):
        headers = ['u_x', 'u_y', 'v_x', 'v_y', 'a_x', 'a_y', '||a||', 'time', 'loadfunc']
        table = {key:[] for key in headers}

        if timelog:
            tlog = open(timelog, 'w')
        else:
            tlog = None

        dt = t_max / float(maxiter)
        t1 = time.clock()
        if timelog: tlog.write('Commencing initializing calculation\n')
        # Initialize local values
        els = self.elements
        nds = self.nodes
        for el in els:
            el.set_matrices()
            el.calc_normal_vectors()
        for nd in nds:
            nd.init_vals(maxiter)
        t2 = time.clock()
        if timelog: tlog.write('Initializing calculation total time = %g\n' % (t2 - t1))

        # Iterate until convergence is reached for all elements
        itercount = 0
        crits = []
        node_data = []
        while True:
            itercount += 1
            crit = 0.
            if timelog: tlog.write('Starting iteration no. %d\n' % itercount)
            t1 = time.clock()
            for i, el in enumerate(els):
                if timelog: tlog.write('Calculating element no. %d\n' % (i + 1))
                crit = max(crit, el.iterate(dt, itercount * dt / t_load, tlog, verbose=verbose))
            for i, el in enumerate(els):
                if timelog: tlog.write('Writing nodal values element no. %d\n' % (i + 1))
                el.write_next_vals()
            t2 = time.clock()
            if timelog: tlog.write('Iteration no. %d complete, total time taken = %g\n' % (itercount, t2 - t1))
            crits.append(crit)
            if node_to_plot and ((itercount % output_freq) == 0):
                node_data.append((np.array(self.nodes[node_to_plot - 1].v_disp),
                                               np.array(self.nodes[node_to_plot - 1].v_velo),
                                               np.array(self.nodes[node_to_plot - 1].v_acce),
                                               ))
                table['||a||'].append(np.linalg.norm(np.hstack([nd.v_acce for nd in nds])))
                table['loadfunc'].append(load_function(itercount * dt / t_load))
                table['time'].append(itercount * dt)
            if crit < tol:
                if verbose: print 'Convergence reached in step no.%d: maximum criterium = %.10f' % (itercount, crit)
                # break
            else:
                if verbose: print 'Step no.%d: maximum criterium = %.10f' % (itercount, crit)
            if itercount == maxiter:
                print 'Maximum step count reached: maximum criterium = %.10f' % crit
                break
        acce_norm = np.linalg.norm(np.hstack([nd.v_acce for nd in nds]))
        print 'Acceleration norm upon calculation end = %.6f' % acce_norm
        if timelog: tlog.close()
        t_vals = np.array(table['time'])
        # plt.plot(t_vals, crits, label='Convergence criterium', linewidth=2.)
        if node_to_plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(t_vals, [n[0][0] for n in node_data], label='u_x', linewidth=2.)
            ax.plot(t_vals, [n[0][1] for n in node_data], label='u_y', linewidth=2.)
            ax.plot(t_vals, [n[1][0] for n in node_data], label='v_x', linewidth=2.)
            ax.plot(t_vals, [n[1][1] for n in node_data], label='v_y', linewidth=2.)
            ax.plot(t_vals, [n[2][0] for n in node_data], label='a_x', linewidth=2.)
            ax.plot(t_vals, [n[2][1] for n in node_data], label='a_y', linewidth=2.)
            table['u_x'] = [n[0][0] for n in node_data]
            table['u_y'] = [n[0][1] for n in node_data]
            table['v_x'] = [n[1][0] for n in node_data]
            table['v_y'] = [n[1][1] for n in node_data]
            table['a_x'] = [n[2][0] for n in node_data]
            table['a_y'] = [n[2][1] for n in node_data]
            table['time'] = t_vals
            ax.legend(loc='lower left')
            if savetextfile:
                f = open(savetextfile, 'w')
                texttable = tbl.tabulate(zip(*[table[col] for col in headers]), headers=headers, tablefmt="orgtbl")
                f.write(texttable)
                f.close()
            if savefile:
                fig.savefig(savefile)
            if show:
                plt.show()

    def solve_kin_damping(self, dt=0.01, tol=0.001, maxiter=50, maxiniter=1000, verbose=False, savefile='', show=True):
        # Initialize local values
        els = self.elements
        nds = self.nodes
        for el in els:
            el.set_matrices()
            el.calc_normal_vectors()
        for nd in nds:
            nd.init_vals(maxiter)

        # Iterate until convergence is reached for all elements
        itercount = 0
        cycles = []
        energies = []
        next_disp, next_velo, next_acce = [[], [], []]
        while True:
            done = False
            diverg = False
            itercount += 1
            if verbose: print 'step %d' % itercount
            initers = 0
            while True:
                initers += 1
                if verbose: print '    inner step %d' % initers
                energy = 0.
                for el in els:
                    iter = el.iterate_kin_damping(dt, itercount / float(maxiter))  # ,verbose
                    energy += iter[0]
                    next_disp.append(iter[1])
                    next_velo.append(iter[2])
                    next_acce.append(iter[3])
                energies.append(energy)

                if initers == maxiniter:
                    print 'Maximum inner step count reached: maximum energy = %.10f, try setting higher inner step count' % energy
                    print '%d outer steps completed' % itercount
                    cycles.append((energies, dt))
                    diverg = True
                    break
                if (len(energies) > 2) and (energies[-1] < energies[-2] > energies[-3]):
                    # Maximum energy reached in previous step, reset velocities
                    for nd in nds:
                        nd.set_velo(np.zeros(2))
                        nd.init_vals(maxiter)
                    if verbose: print '    Energy: ', energy
                    cycles.append((energies, dt))
                    energies = []
                    dt *= 0.85
                    next_disp, next_velo, next_acce = [[], [], []]
                    if energy < tol: done = True
                    break
                else:
                    # Maximum energy not reached, continue iterating
                    for i, el in enumerate(els):
                        el.set_node_vals(next_disp[i], next_velo[i], next_acce[i])
                    next_disp, next_velo, next_acce = [[], [], []]
                    if verbose: print '    Energy: ', energy
            if done:
                print 'Convergence reached in %d steps, E=%f' % (itercount, energy)
                break
            if diverg or itercount > maxiter:
                break
        if show or savefile:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            e_lst = [np.array(c[0]) for c in cycles]
            e = np.hstack(e_lst)
            t_lst = []
            t0 = 0.
            for c in cycles:
                t_lst.append(np.linspace(t0 + c[1], t0 + len(c[0]) * c[1], len(c[0])))
                t0 += t_lst[-1][-1]
            t = np.hstack(t_lst)
            ax.plot(t, e, label='Kinetic energy', linewidth=2.)
            ax.legend(loc='lower left')
            if savefile:
                fig.savefig(savefile)
            if show:
                plt.show()

class DomainFEM(Domain):
    def solve(self, verbose=False):
        # Initialize local values
        els = self.elements
        nds = self.nodes
        code_count = 0
        for nd in nds:
            code_count = nd.set_codes(code_count, verbose)
        if verbose: print 'Number of unknowns: %d' % code_count

        f_glob = np.zeros(code_count)  # Load vector
        for nd in nds:
            for i, c in enumerate(nd.v_code):
                if c != 0: f_glob[c - 1] += nd.F_ext[i]

        K_glob = np.zeros((code_count, code_count))  # Global stiffness
        for el in els:
            el.set_matrices()
            el.set_codes()
            cds = [(i, el.v_code[i]) for i in range(len(el.v_code)) if el.v_code[i] != 0]
            K_loc = el.K
            for ii, ij in it.product(cds, cds):
                K_glob[ii[1] - 1, ij[1] - 1] += K_loc[ii[0], ij[0]]

        K_inv = np.linalg.inv(K_glob)
        u_res = np.dot(K_inv, f_glob)
        if verbose:
            print 'Global stiffness matrix:', K_glob
            print 'Inverse stiffness matrix:', K_inv
            print "Global load vector", f_glob
            print "Resulting displacement vector", u_res
        for nd in nds:
            for i, c in enumerate(nd.v_code):
                nd.v_disp[i] = (c != 0) * u_res[c - 1]

class DomainCDEMStatic(Domain):
    def __init__(self, c_n, c_s, elements=[], nodes=[], stiffness_scaling=False):
        self.elements = elements
        for el in self.elements:
            el.domain = weakref.ref(self)
        self.nodes = nodes
        self.stiffness_scaling = stiffness_scaling  # Contact forces are multiplied by edge length
        self.set_contacts(c_n, c_s)

    def set_contacts(self, c_n, c_s):
        self.C = np.array([[c_n, 0.],
                           [0., c_s]])  # Contact stiffness matrix

    def solve(self, verbose=False, timelog=None):
        # Initialize local values
        els = self.elements
        nds = self.nodes

        t1 = time.clock()
        code_count = 0
        for nd in nds:
            code_count = nd.set_codes(code_count, verbose)
        t2 = time.clock()
        if verbose: print 'Number of unknowns: %d' % code_count

        K_glob = np.zeros((code_count, code_count))  # Global stiffness
        for el in els:
            el.set_matrices()
            el.set_codes()
            el.calc_normal_vectors()
            cds = [(i, el.v_code[i]) for i in range(len(el.v_code)) if el.v_code[i] != 0]
            K_loc = el.K
            for ii, ij in it.product(cds, cds):
                K_glob[ii[1] - 1, ij[1] - 1] += K_loc[ii[0], ij[0]]
        t3 = time.clock()

        f_glob = np.zeros(code_count)  # Load vector
        for nd in nds:
            for i, c in enumerate(nd.v_code):
                if c != 0: f_glob[c - 1] += (c != 0) * nd.F_ext[i]

            # Contact components of the global stiffness matrix
            for nbr, n_vect, l_edge in zip(nd.neighbors, nd.n_vects, nd.l_edges):
                if nbr:
                    T = np.array([[n_vect[0], n_vect[1]],
                                  [-n_vect[1], n_vect[0]]])  # Transformation matrix
                    C = np.dot(np.dot(np.transpose(T), self.C), T)
                    if self.stiffness_scaling: C *= l_edge / 2.
                    c1 = [(i, nd.v_code[i]) for i in range(len(nd.v_code)) if nd.v_code[i] != 0]
                    nbrnd = nds[nbr - 1]
                    c2 = [(i, nbrnd.v_code[i]) for i in range(len(nbrnd.v_code)) if nbrnd.v_code[i] != 0]
                    for i, j in it.product(c1, c1 + c2):
                        # du = u_i - u_0
                        if i[1] == j[1]:
                            K_glob[i[1] - 1, j[1] - 1] += C[i[0], j[0]]
                        else:
                            K_glob[i[1] - 1, j[1] - 1] -= C[i[0], j[0]]
        t4 = time.clock()

        K_inv = np.linalg.inv(K_glob)
        t5 = time.clock()
        u_res = np.dot(K_inv, f_glob)
        t6 = time.clock()
        for nd in nds:
            for i, c in enumerate(nd.v_code):
                nd.v_disp[i] = (c != 0) * u_res[c - 1]
        t7 = time.clock()

        if timelog:
            tlog = open(timelog, 'w')
            tlog.write('    Code assigning: %g\n' % (t2 - t1))
            tlog.write('    Element global stiffness: %g\n' % (t3 - t2))
            tlog.write('    load vector + contact global stiffness: %g\n' % (t4 - t3))
            tlog.write('    Inverting global stiffness matrix: %g\n' % (t5 - t4))
            tlog.write('    Calculating result: %g\n' % (t6 - t5))
            tlog.write('    Distributing result: %g\n' % (t7 - t6))
            tlog.close()

        if verbose:
            print 'Global stiffness matrix:', K_glob
            print 'Inverse stiffness matrix:', K_inv
            print "Global load vector", f_glob
            print "Resulting displacement vector", u_res

if __name__ == '__main__':
    pass
