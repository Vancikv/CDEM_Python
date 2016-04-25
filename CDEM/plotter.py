'''
Created on 14. 1. 2016

@author: Kapsak
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import weakref

class PlotNode(object):
    def __init__(self, x, y, v_disp, v_velo, v_acce):
        self.x = x
        self.y = y
        self.v_disp = np.array([v_disp])
        self.v_velo = np.array([v_velo])
        self.v_acce = np.array([v_acce])

class PlotElement(object):
    def __init__(self, nodes, ref):
        self.nodes = nodes
        self.ref = ref

    def plot(self, n_step, magnitude, linestyle='-', color='blue', linewidth=3.):
        ref = self.ref()
        nds = self.nodes + [self.nodes[0]]
        xy = np.transpose(np.array([[ref.nodes[i - 1].x, ref.nodes[i - 1].y] for i in nds]))
        u = np.transpose(np.array([ref.nodes[i - 1].v_disp[n_step] for i in nds]))

        ref.ax.plot(xy[0, :] + magnitude * u[0, :], xy[1, :] + magnitude * u[1, :], linestyle=linestyle, color=color, linewidth=linewidth)


class ResultPlotter(object):
    def __init__(self):
        self.nodes = []
        self.elements = []
        self.fig = plt.figure()
        self.ax = self.fig.add_axes([.2, .2, .7, .7])

    def load_data(self, path, prefix):
        lstdir = os.listdir(path)
        files = []
        for file in lstdir:
            if os.path.split(file)[1].find(prefix) == 0:
                files.append(os.path.join(path, file))
        time_dict = {}
        for file in files:
            fl = open(file, 'r')
            ln = fl.readline()
            time_dict[float(ln.split(' ')[2])] = file
            fl.close()
        timeline = sorted(time_dict.keys())
        self.timeline = timeline

        fl = open(time_dict[timeline[0]], 'r')
        ln = fl.readline().split(' ')
        nnodes = int(ln[0])
        nelems = int(ln[1])
        loadfunc = [float(ln[3])]
        for i in range(nnodes):
            ln = fl.readline().split(' ')
            x, y = float(ln[2]), float(ln[3])
            v_disp = [float(ln[4]), float(ln[5])]
            v_velo = [float(ln[6]), float(ln[7])]
            v_acce = [float(ln[8]), float(ln[9])]
            self.nodes.append(PlotNode(x, y, v_disp, v_velo, v_acce))
        for i in range(nelems):
            ln = fl.readline().split(' ')
            self.elements.append(PlotElement(map(int, ln[3:]), weakref.ref(self)))
        fl.close()

        for key in timeline[1:]:
            fl = open(time_dict[key], 'r')
            ln = fl.readline().split(' ')
            loadfunc.append(float(ln[3]))
            for i in range(nnodes):
                ln = fl.readline().split(' ')
                v_disp = [float(ln[4]), float(ln[5])]
                v_velo = [float(ln[6]), float(ln[7])]
                v_acce = [float(ln[8]), float(ln[9])]
                self.nodes[i].v_disp = np.vstack((self.nodes[i].v_disp, v_disp))
                self.nodes[i].v_velo = np.vstack((self.nodes[i].v_velo, v_velo))
                self.nodes[i].v_acce = np.vstack((self.nodes[i].v_acce, v_acce))
            fl.close()
        self.loadfunc = loadfunc

    def save_plot(self, path):
        self.fig.savefig(path)
        self.ax.cla()
        pass

    def plot_step(self, n_step, magnitude=1., linestyle='-', color='blue', linewidth=3.):
        for e in self.elements:
            e.plot(n_step, magnitude, linestyle=linestyle, color=color, linewidth=linewidth)
        self.ax.axis('equal')

    def plot_loadfunc(self, t1, t2, label, xlbl='', ylbl=''):
        if t1 > t2:
            print "Error in plot_graph, wrong time boundaries (t1 > t2)."
        if t1 < 0.:
            print "Error in plot_graph, wrong time boundaries (t1 < 0.)."
        tl = self.timeline

        # Obtain step boundaries from time boundaries.
        s1, s2 = 0, len(tl)
        for i, t in enumerate(tl):
            if t < t1: s1 += 1
            if t > t2:
                s2 = i
                break

        plot_timeline = tl[s1:s2]
        plot_vals = self.loadfunc[s1:s2]

        self.ax.plot(plot_timeline, plot_vals, label=label)
        self.ax.set_ylim([-1, 1])
        self.ax.legend()
        self.ax.set_xlabel(xlbl, fontsize=24)
        self.ax.set_ylabel(ylbl, fontsize=24)

    def plot_graph(self, t1, t2, val_name, val_id, node, label, xlbl='', ylbl=''):
        if t1 > t2:
            print "Error in plot_graph, wrong time boundaries (t1 > t2)."
        if t1 < 0.:
            print "Error in plot_graph, wrong time boundaries (t1 < 0.)."
        tl = self.timeline

        # Obtain step boundaries from time boundaries.
        s1, s2 = 0, len(tl)
        for i, t in enumerate(tl):
            if t < t1: s1 += 1
            if t > t2:
                s2 = i
                break

        plot_timeline = tl[s1:s2]
        plot_node = self.nodes[node - 1]
        plot_vals = getattr(plot_node, val_name)[s1:s2, val_id]

        self.ax.plot(plot_timeline, plot_vals, label=label)
        self.ax.legend()
        self.ax.set_xlabel(xlbl, fontsize=24)
        self.ax.set_ylabel(ylbl, fontsize=24)
        pass

a = ResultPlotter()

# a.load_data("C:/Users/Werner/Dropbox/VladimirVancik/CUDA/jobs/3pb", "3pbo")
# a.plot_graph(0., 30.0, 'v_disp', 1, 59, 'v_disp', xlbl=r'$t$', ylbl=r'$u_y$')
# a.save_plot("C:/Users/Werner/Dropbox/VladimirVancik/CUDA/jobs/3pb/3pbgraph.png")
# a.plot_step(119, 1000000)
# a.save_plot("C:/Users/Werner/Dropbox/VladimirVancik/CUDA/jobs/3pb/3pbstep100000pic.png")
# a.load_data("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb", "3pbout")
# a.load_data("C:/Users/Kapsak/Dropbox/VladimirVancik/CUDA/jobs/3pb", "3pbo")
# a.plot_graph(0., 2.0, 'v_disp', 1, 59, 'v_disp', xlbl=r'$t$', ylbl=r'$u_y$')
# a.save_plot("C:/Users/Kapsak/Dropbox/VladimirVancik/CUDA/jobs/3pb/3pbgraph.png")
# a.plot_step(99, 100000)
# a.save_plot("C:/Users/Kapsak/Dropbox/VladimirVancik/CUDA/jobs/3pb/3pbstep100000pic.png")
# a.load_data("C:/Users/Kapsak/Dropbox/CDEM/Cpp_data/3pb", "3pbout")
# a.plot_loadfunc(0., 30.0, 'Load function', xlbl=r'$t [s]$', ylbl=r'$lf [-]$')
# a.save_plot("C:/Users/Kapsak/Dropbox/CDEM/Cpp_data/3pb/3pb_loadfunc.png")
# a.plot_graph(0., 30.0, 'v_disp', 1, 59, 'v_disp', xlbl=r'$t$', ylbl=r'$u_y$')
# a.save_plot("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb/3pbgraph_new.png")
# a.plot_step(119, 100)
# a.save_plot("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb/3pbendpic_new.png")
a.load_data("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb_long", "3pbout")
a.plot_loadfunc(0.0, 0.03, 'Load function', xlbl=r'$t [s]$', ylbl=r'$lf [-]$')
a.save_plot("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb_long/3pb_loadfunc.pdf")
a.plot_graph(0.0, 0.03, 'v_disp', 1, 119, 'Center deflection', xlbl=r'$t$', ylbl=r'$u_y$')
a.save_plot("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb_long/3pbgraph.pdf")
a.plot_step(120, 10)
a.save_plot("C:/Users/Werner/Dropbox/CDEM/Cpp_data/3pb_long/3pbstep10000pic.pdf")

