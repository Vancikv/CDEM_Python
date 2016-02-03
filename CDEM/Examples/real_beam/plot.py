'''
Created on Dec 14, 2015

@author: vancik
'''

import numpy as np
import os
import matplotlib.pyplot as plt

def get_outfile(folder_path_list, file_name):
    '''Returns a file in the specified folder using the home
    directory as root.
    '''
    HOME_DIR = os.path.expanduser("~")
    out_dir = os.path.join(HOME_DIR, *folder_path_list)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    outfile = os.path.join(out_dir, file_name)
    return outfile

f = plt.figure(figsize=(10, 7))
ax = f.add_axes([0.12, 0.12, 0.84, 0.78])
dataid = 0
mini = 1e10
maxi = -1e10
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbe_004.txt'))
ax.plot(np.linspace(0., 0.004, 40), data[:, dataid], label=r'OOFEM EX, $dt=5e-6$', linestyle='-', color='black')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbe_008.txt'))
ax.plot(np.linspace(0., 0.008, 80), data[:, dataid], label=None, linestyle='-', color='black')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbi_004.txt'))
ax.plot(np.linspace(0., 0.004, 40), data[:, dataid], label=r'OOFEM IM, $dt=1e-4$', linestyle='-', color='blue')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbi_008.txt'))
ax.plot(np.linspace(0., 0.008, 80), data[:, dataid], label=None, linestyle='-', color='blue')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=160000.0.txt'), delimiter='|', skiprows=3, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 1], label=r'Py EX, $dt=2.5e-8$', linestyle='-.', color='red')
mini = min(mini, np.min(data[:, 1]))
maxi = max(maxi, np.max(data[:, 1]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=320000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 1], label=r'Py EX, $dt=1.25e-8$', linestyle='--', color='red')
mini = min(mini, np.min(data[:, 1]))
maxi = max(maxi, np.max(data[:, 1]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=640000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 1], label=r'Py EX, $dt=0.625e-8$', linestyle='-', color='red')
mini = min(mini, np.min(data[:, 1]))
maxi = max(maxi, np.max(data[:, 1]))
# data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=1280000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
# ax.plot(data[:, 7], data[:, 1], label=r'Py EX, $dt=0.3125e-8$', linestyle='-', color='red')
# mini = min(mini, np.min(data[:, 1]))
# maxi = max(maxi, np.max(data[:, 1]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_008'], 'textgraph_niter=320000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 1], label=None, linestyle='-.', color='red')
mini = min(mini, np.min(data[:, 1]))
maxi = max(maxi, np.max(data[:, 1]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_008'], 'textgraph_niter=640000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 1], label=None, linestyle='--', color='red')
mini = min(mini, np.min(data[:, 1]))
maxi = max(maxi, np.max(data[:, 1]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_008'], 'textgraph_niter=1280000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 1], label=None, linestyle='-', color='red')
mini = min(mini, np.min(data[:, 1]))
maxi = max(maxi, np.max(data[:, 1]))
ax.set_ylabel(r'$u_y [m]$', fontsize=18)
ax.set_xlabel(r'$t [s]$', fontsize=18)
ax.set_title(r'Vertical displacement of center node, load application times of $0.004s$ and $0.008s$', y=1.05, fontsize=16)
ax.set_ylim((mini - 0.15 * abs(mini), 1.15 * maxi))
ax.legend(loc='lower left', fontsize=14)
zed = [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
zed = [tick.label.set_fontsize(14) for tick in ax.xaxis.get_major_ticks()]
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
f.savefig(get_outfile(['Dropbox', 'CDEM', '3point_bending_real'], 'uy_gradual_t004x008.png'))
plt.cla()

####################
### ACCELERATION ###
####################

dataid = 1
mini = 1e10
maxi = -1e10
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbe_004.txt'))
ax.plot(np.linspace(0., 0.004, 40), data[:, dataid], label=r'OOFEM EX, $dt=5e-6$', linestyle='-', color='black')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbe_008.txt'))
ax.plot(np.linspace(0., 0.008, 80), data[:, dataid], label=None, linestyle='-', color='black')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbi_004.txt'))
ax.plot(np.linspace(0., 0.004, 40), data[:, dataid], label=r'OOFEM IM, $dt=1e-4$', linestyle='-', color='blue')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', 'OOFEM_input', 'real_beam', 'porovnani'], '3pbi_008.txt'))
ax.plot(np.linspace(0., 0.008, 80), data[:, dataid], label=None, linestyle='-', color='blue')
mini = min(mini, np.min(data[:, dataid]))
maxi = max(maxi, np.max(data[:, dataid]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=160000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 5], label=r'Py EX, $dt=2.5e-8$', linestyle='-.', color='red')
mini = min(mini, np.min(data[:, 5]))
maxi = max(maxi, np.max(data[:, 5]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=320000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 5], label=r'Py EX, $dt=1.25e-8$', linestyle='--', color='red')
mini = min(mini, np.min(data[:, 5]))
maxi = max(maxi, np.max(data[:, 5]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=640000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 5], label=r'Py EX, $dt=0.625e-8$', linestyle='-', color='red')
mini = min(mini, np.min(data[:, 5]))
maxi = max(maxi, np.max(data[:, 5]))
# data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_004'], 'textgraph_niter=1280000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
# ax.plot(data[:, 7], data[:, 5], label=r'Py EX, $dt=0.3125e-8$', linestyle='-', color='red')
# mini = min(mini, np.min(data[:, 5]))
# maxi = max(maxi, np.max(data[:, 5]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_008'], 'textgraph_niter=320000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 5], label=None, linestyle='-.', color='red')
mini = min(mini, np.min(data[:, 5]))
maxi = max(maxi, np.max(data[:, 5]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_008'], 'textgraph_niter=640000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 5], label=None, linestyle='--', color='red')
mini = min(mini, np.min(data[:, 5]))
maxi = max(maxi, np.max(data[:, 5]))
data = np.loadtxt(get_outfile(['Dropbox', 'CDEM', '3point_bending_real', 'totaltime_008'], 'textgraph_niter=1280000.0.txt'), delimiter='|', skiprows=2, usecols=range(1, 10))
ax.plot(data[:, 7], data[:, 5], label=None, linestyle='-', color='red')
mini = min(mini, np.min(data[:, 5]))
maxi = max(maxi, np.max(data[:, 5]))
ax.set_ylabel(r'$a_y [ms^{-2}]$', fontsize=18)
ax.set_xlabel(r'$t [s]$', fontsize=18)
ax.set_title(r'Acceleration of center node, load application times of $0.004s$ and $0.008s$', y=1.05, fontsize=16)
ax.set_ylim((mini - 0.15 * abs(mini), 1.15 * maxi))
ax.legend(loc='lower left', fontsize=14)
zed = [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
zed = [tick.label.set_fontsize(14) for tick in ax.xaxis.get_major_ticks()]
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
f.savefig(get_outfile(['Dropbox', 'CDEM', '3point_bending_real'], 'ay_gradual_t004x008.png'))
plt.cla()

