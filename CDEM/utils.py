'''
Created on 10. 11. 2015

@author: Kapsak
'''

from basic import *
import numpy as np
import weakref
import copy
from domains import *
from nodes import *
from elements import *

type_dict = {'CDEMDR':DomainCDEMDR, 'CDEMStatic':DomainCDEMStatic}

def get_2D_quadrangle_domain(ni, nj, lx, ly, E=0., nu=0., density=0.,
                             thickness=1., alfaC=1., stiffdef=None,
                             c_n=1000., c_s=1000.,
                             load=[], supports=[],
                             domtype='CDEMDR', spring_stiffness_factor=0.):
    '''
    Return a domain with an orthogonal mesh.
    ni = elements in y direction, nj = elements in x direction
    lx = length of specimen, ly = height of specimen
    dx = element size in x direction, dy = element size in y direction

    load and support format:
    [ ((xmin, xmax), (ymin, ymax), val_lst) ]
    '''
    dx = lx / float(nj)
    dy = ly / float(ni)
    nodes = []
    elements = []
    dom = type_dict[domtype](c_n=c_n, c_s=c_s)
    for j in range(nj):
        for i in range(ni):
            lsx = j * dx * np.ones(4) + np.array([0., dx, dx, 0.])
            lsy = i * dy * np.ones(4) + np.array([0., 0., dy, dy])
            for x, y in zip(lsx, lsy):
                nodes.append(Node(x=x, y=y, neighbors=[0, 0]))

                # Check if load or support should be assigned
                for bx, by, val in load:
                    if all((x > bx[0], x < bx[1], y > by[0], y < by[1])):
                        nodes[-1].F_ext = val
                for bx, by, val in supports:
                    if all((x > bx[0], x < bx[1], y > by[0], y < by[1])):
                        nodes[-1].supports = val
            if i > 0:
                nodes[-4].neighbors[1] = j * ni * 4 + 4 * (i - 1) + 4
                nodes[-3].neighbors[0] = j * ni * 4 + 4 * (i - 1) + 3
            if i < ni - 1:
                nodes[-2].neighbors[1] = j * ni * 4 + 4 * (i + 1) + 2
                nodes[-1].neighbors[0] = j * ni * 4 + 4 * (i + 1) + 1
            if j > 0:
                nodes[-4].neighbors[0] = (j - 1) * ni * 4 + 4 * i + 2
                nodes[-1].neighbors[1] = (j - 1) * ni * 4 + 4 * i + 3
            if j < nj - 1:
                nodes[-3].neighbors[1] = (j + 1) * ni * 4 + 4 * i + 1
                nodes[-2].neighbors[0] = (j + 1) * ni * 4 + 4 * i + 4
            idx = j * ni * 4 + 4 * i
            nodeids = [idx + 1, idx + 2, idx + 3, idx + 4]
            elements.append(ElemQuadrangleLin(E=E, nu=nu, density=density, thickness=thickness,
                                              domain=weakref.ref(dom), alfaC=alfaC, stiffdef=stiffdef, nodes=nodeids))
    dom.nodes = nodes
    dom.elements = elements
    if spring_stiffness_factor:
        max_k = 0.
        for el in elements:
            if not hasattr(el, 'K'):
                el.set_matrices()
            max_k = max(max_k, np.max(el.K))
        c_n = max_k * spring_stiffness_factor
        c_s = max_k * spring_stiffness_factor
        dom.set_contacts(c_n, c_s)
    return dom

def get_2D_quadrangle_FEM_domain(ni, nj, lx, ly, E=0., nu=0., density=0.,
                             thickness=1., stiffdef=None,
                             load=[], supports=[]):
    '''
    Return a domain with an orthogonal mesh.
    ni = elements in y direction, nj = elements in x direction
    lx = length of specimen, ly = height of specimen
    dx = element size in x direction, dy = element size in y direction

    load and support format:
    [ ((xmin, xmax), (ymin, ymax), val_lst) ]
    '''
    dx = lx / float(nj)
    dy = ly / float(ni)
    nodes = []
    elements = []
    dom = DomainFEM()
    for j in range(nj + 1):
        for i in range(ni + 1):
            nodes.append(Node(x=j * dx, y=i * dy, neighbors=[0, 0]))
            # Check if load or support should be assigned
            x, y = nodes[-1].x, nodes[-1].y
            for bx, by, val in load:
                if all((x > bx[0], x < bx[1], y > by[0], y < by[1])):
                    nodes[-1].F_ext = val
            for bx, by, val in supports:
                if all((x > bx[0], x < bx[1], y > by[0], y < by[1])):
                    nodes[-1].supports = val
            if i < ni and j < nj:
                idx = j * (ni + 1) + i
                nodeids = [idx + 1, idx + 2 + ni, idx + 3 + ni, idx + 2]
                elements.append(ElemQuadrangleLin(E=E, nu=nu, density=density, thickness=thickness,
                                                  domain=weakref.ref(dom), stiffdef=stiffdef, nodes=nodeids))
    dom.nodes = nodes
    dom.elements = elements
    return dom

def get_stencil_mesh(stencil, dx, dy, x0, y0, n0):
    if stencil == '0':
        x = [x0 + dxi for dxi in [0., dx, dx, 0.]]
        y = [y0 + dyi for dyi in [0., 0., dy, dy]]
        nodes = [Node(x=xi, y=yi, neighbors=[0, 0]) for xi, yi in zip(x, y)]
        elements = [ElemQuadrangleLin(nodes=[n0 + i + 1 for i in range(4)])]
    if stencil == 'f':
        x = [x0 + dxi for dxi in [0., dx, dx, 0., dx, dx, 0., 0., dx, dx]]
        y = [y0 + dyi for dyi in [0., 0., dy / 3., 0., dy / 3., 2 * dy / 3., dy, dy, 2 * dy / 3., dy]]
        nodes = [Node(x=xi, y=yi, neighbors=[0, 0]) for xi, yi in zip(x, y)]

        # Internal node connectivity
        nodes[0].neighbors = [n0 + 4, 0]
        nodes[2].neighbors = [0, n0 + 5]
        nodes[3].neighbors = [0, n0 + 1]
        nodes[4].neighbors = [n0 + 3, 0]
        nodes[5].neighbors = [0, n0 + 9]
        nodes[6].neighbors = [n0 + 8, 0]
        nodes[7].neighbors = [0, n0 + 7]
        nodes[8].neighbors = [n0 + 6, 0]

        elements = [ElemTriangleLin(nodes=[n0 + i + 1 for i in range(3)]),
                    ElemQuadrangleLin(nodes=[n0 + i + 1 for i in range(3, 7)]),
                    ElemTriangleLin(nodes=[n0 + i + 1 for i in range(7, 10)])]
    if stencil == 'c':
        x = [x0 + dxi for dxi in [0., dx, 0., 0., dx, dx, 0., 0., dx, 0.]]
        y = [y0 + dyi for dyi in [0., 0., dy / 3., dy / 3., 0., dy, 2 * dy / 3., 2 * dy / 3., dy, dy]]
        nodes = [Node(x=xi, y=yi, neighbors=[0, 0]) for xi, yi in zip(x, y)]

        # Internal node connectivity
        nodes[1].neighbors = [0, n0 + 5]
        nodes[2].neighbors = [n0 + 4, 0]
        nodes[3].neighbors = [0, n0 + 3]
        nodes[4].neighbors = [n0 + 2, 0]
        nodes[5].neighbors = [0, n0 + 9]
        nodes[6].neighbors = [n0 + 8, 0]
        nodes[7].neighbors = [0, n0 + 7]
        nodes[8].neighbors = [n0 + 6, 0]

        elements = [ElemTriangleLin(nodes=[n0 + i + 1 for i in range(3)]),
                    ElemQuadrangleLin(nodes=[n0 + i + 1 for i in range(3, 7)]),
                    ElemTriangleLin(nodes=[n0 + i + 1 for i in range(7, 10)])]
    return (nodes, elements)

def get_uneven_mesh(lx, ly, nj, ni0, stencil_arr,
                             E=0., nu=0., density=0.,
                             thickness=1., alfaC=1., stiffdef=None,
                             c_n=100000., c_s=100000.,
                             load=[], supports=[],
                             domtype='CDEMDR', spring_stiffness_factor=0.):
    '''
    Return a CDEM domain with a rectangular mesh that can grow finer or coarser
    along the width.

    nj - number of vertical element groups (length of the stencil_arr)
    ni0 - number of element stencils in the first vertical group
    stencil_arr - each entry specifies the stencil type for one vertical
                  element group, left to right, length has to be nj
                  0 - same element size going right
                  f - 3x finer mesh going right
                  c - 3x coarser mesh going right
    '''
    number_of_nodes = {'0':4, 'f':10, 'c':10}
    number_of_elements = {'0':1, 'f':3, 'c':3}
    arg_ni = {'0':lambda ni: ni, 'f':lambda ni: ni * 3, 'c':lambda ni: ni}
    get_ni = {'0':lambda ni: ni, 'f':lambda ni: ni, 'c':lambda ni: ni / 3}
    arg_dx = {'0':lambda dx: dx, 'f':lambda dx: dx / 2., 'c':lambda dx: dx * 3. / 2.}
    get_dx = {'0':lambda dx: dx, 'f':lambda dx: 2. / 3.*dx, 'c':lambda dx: dx * 2.}
    dx, dy, x0 = [np.zeros(nj) for i in range(3)]
    ni, n0, e0 = [np.zeros(nj, dtype=int) for i in range(3)]
    ni[0] = get_ni[stencil_arr[0]](ni0)  # number of stencil instances in vertical group
    dy[0] = ly / float(ni[0])  # y length of vertical group
    dx[0] = 1.  # x length of vertical group (unit at first)
    n0[0] = 0  # id of first node of vertical group
    x0[0] = 0.  # x coordinate of first node of vertical group
    e0[0] = 0  # id of first element of vertical group
    for i in range(nj - 1):
        ni[i + 1] = get_ni[stencil_arr[i + 1]](arg_ni[stencil_arr[i]](ni[i]))
        dx[i + 1] = get_dx[stencil_arr[i + 1]](arg_dx[stencil_arr[i]](dx[i]))
        dy[i + 1] = ly / float(ni[i + 1])
        n0[i + 1] = n0[i] + ni[i] * number_of_nodes[stencil_arr[i]]
        x0[i + 1] = x0[i] + dx[i]
        e0[i + 1] = e0[i] + ni[i] * number_of_elements[stencil_arr[i]]
    lx_unit = lx / sum(dx)
    dx *= lx_unit  # assign real length
    x0 *= lx_unit

    # Create nodes and elements
    nodes = []
    elements = []
    for j in range(nj):
        for i in range(ni[j]):
            stencil_mesh = get_stencil_mesh(stencil_arr[j], dx[j], dy[j], x0[j], i * dy[j], n0[j] + i * number_of_nodes[stencil_arr[j]])
            nodes += stencil_mesh[0]
            elements += stencil_mesh[1]

    # Make connections
    for j in range(nj):
        for i in range(ni[j]):
            # Internal connections are already present
            contact = ''.join(stencil_arr[j:j + 2])  # Store the type of contact
            ns1 = number_of_nodes[stencil_arr[j]]
            nbl1 = n0[j] + i * ns1  # bottom left corner of bottom left stencil

            # j-external connections (right stencil edge)
            if nj - 1 > j:
                ns2 = number_of_nodes[stencil_arr[j + 1]]
                if contact == 'ff':
                    nbl2 = n0[j + 1] + 3 * i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 1].neighbors[1] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[0] = nbl1 + 2
                    nodes[nbl1 + 2].neighbors[0] = nbl2 + 7
                    nodes[nbl2 + 6].neighbors[1] = nbl1 + 3

                    nbl2 += ns2
                    nodes[nbl1 + 4].neighbors[1] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[0] = nbl1 + 5
                    nodes[nbl1 + 5].neighbors[0] = nbl2 + 7
                    nodes[nbl2 + 6].neighbors[1] = nbl1 + 6

                    nbl2 += ns2
                    nodes[nbl1 + 8].neighbors[1] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[0] = nbl1 + 9
                    nodes[nbl1 + 9].neighbors[0] = nbl2 + 7
                    nodes[nbl2 + 6].neighbors[1] = nbl1 + 10
                elif contact == 'fc':
                    nbl2 = n0[j + 1] + i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 1].neighbors[1] = nbl2 + 1
                    nodes[nbl2].neighbors[0] = nbl1 + 2
                    nodes[nbl1 + 2].neighbors[0] = nbl2 + 3
                    nodes[nbl2 + 2].neighbors[1] = nbl1 + 3

                    nodes[nbl1 + 4].neighbors[1] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[0] = nbl1 + 5
                    nodes[nbl1 + 5].neighbors[0] = nbl2 + 7
                    nodes[nbl2 + 6].neighbors[1] = nbl1 + 6

                    nodes[nbl1 + 8].neighbors[1] = nbl2 + 8
                    nodes[nbl2 + 7].neighbors[0] = nbl1 + 9
                    nodes[nbl1 + 9].neighbors[0] = nbl2 + 10
                    nodes[nbl2 + 9].neighbors[1] = nbl1 + 10
                elif contact == 'f0':
                    nbl2 = n0[j + 1] + 3 * i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 1].neighbors[1] = nbl2 + 1
                    nodes[nbl2].neighbors[0] = nbl1 + 2
                    nodes[nbl1 + 2].neighbors[0] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[1] = nbl1 + 3

                    nbl2 += ns2
                    nodes[nbl1 + 4].neighbors[1] = nbl2 + 1
                    nodes[nbl2].neighbors[0] = nbl1 + 5
                    nodes[nbl1 + 5].neighbors[0] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[1] = nbl1 + 6

                    nbl2 += ns2
                    nodes[nbl1 + 8].neighbors[1] = nbl2 + 1
                    nodes[nbl2].neighbors[0] = nbl1 + 9
                    nodes[nbl1 + 9].neighbors[0] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[1] = nbl1 + 10
                elif contact == 'cf':
                    nbl2 = n0[j + 1] + i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 4].neighbors[1] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[0] = nbl1 + 5
                    nodes[nbl1 + 5].neighbors[0] = nbl2 + 7
                    nodes[nbl2 + 6].neighbors[1] = nbl1 + 6
                elif contact == 'cc':
                    nbl2 = n0[j + 1] + (i // 3) * ns2  # bottom left corner of bottom right stencil
                    if i % 3 == 0:
                        nodes[nbl1 + 4].neighbors[1] = nbl2 + 1
                        nodes[nbl2].neighbors[0] = nbl1 + 5
                        nodes[nbl1 + 5].neighbors[0] = nbl2 + 3
                        nodes[nbl2 + 2].neighbors[1] = nbl1 + 6
                    if i % 3 == 1:
                        nodes[nbl1 + 4].neighbors[1] = nbl2 + 4
                        nodes[nbl2 + 3].neighbors[0] = nbl1 + 5
                        nodes[nbl1 + 5].neighbors[0] = nbl2 + 7
                        nodes[nbl2 + 6].neighbors[1] = nbl1 + 6
                    if i % 3 == 2:
                        nodes[nbl1 + 4].neighbors[1] = nbl2 + 8
                        nodes[nbl2 + 7].neighbors[0] = nbl1 + 5
                        nodes[nbl1 + 5].neighbors[0] = nbl2 + 10
                        nodes[nbl2 + 9].neighbors[1] = nbl1 + 6
                elif contact == 'c0':
                    nbl2 = n0[j + 1] + i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 4].neighbors[1] = nbl2 + 1
                    nodes[nbl2].neighbors[0] = nbl1 + 5
                    nodes[nbl1 + 5].neighbors[0] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[1] = nbl1 + 6
                elif contact == '0f':
                    nbl2 = n0[j + 1] + i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 1].neighbors[1] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[0] = nbl1 + 2
                    nodes[nbl1 + 2].neighbors[0] = nbl2 + 7
                    nodes[nbl2 + 6].neighbors[1] = nbl1 + 3
                elif contact == '0c':
                    nbl2 = n0[j + 1] + (i // 3) * ns2  # bottom left corner of bottom right stencil
                    if i % 3 == 0:
                        nodes[nbl1 + 1].neighbors[1] = nbl2 + 1
                        nodes[nbl2].neighbors[0] = nbl1 + 2
                        nodes[nbl1 + 2].neighbors[0] = nbl2 + 3
                        nodes[nbl2 + 2].neighbors[1] = nbl1 + 3
                    if i % 3 == 1:
                        nodes[nbl1 + 1].neighbors[1] = nbl2 + 4
                        nodes[nbl2 + 3].neighbors[0] = nbl1 + 2
                        nodes[nbl1 + 2].neighbors[0] = nbl2 + 7
                        nodes[nbl2 + 6].neighbors[1] = nbl1 + 3
                    if i % 3 == 2:
                        nodes[nbl1 + 1].neighbors[1] = nbl2 + 8
                        nodes[nbl2 + 7].neighbors[0] = nbl1 + 2
                        nodes[nbl1 + 2].neighbors[0] = nbl2 + 10
                        nodes[nbl2 + 9].neighbors[1] = nbl1 + 3
                elif contact == '00':
                    nbl2 = n0[j + 1] + i * ns2  # bottom left corner of bottom right stencil
                    nodes[nbl1 + 1].neighbors[1] = nbl2 + 1
                    nodes[nbl2].neighbors[0] = nbl1 + 2
                    nodes[nbl1 + 2].neighbors[0] = nbl2 + 4
                    nodes[nbl2 + 3].neighbors[1] = nbl1 + 3
            # i-external connections (top stencil edge)
            if i < ni[j] - 1:
                # Get node indexes
                tl1 = nbl1 + ns1 - 1
                tr1 = nbl1 + ns1 - 2
                bl2 = nbl1 + ns1
                br2 = nbl1 + ns1 + 1
                if stencil_arr[j] == 'f':
                    tl1 = nbl1 + ns1 - 3
                    tr1 = nbl1 + ns1 - 1
                nodes[tl1].neighbors[0] = bl2 + 1
                nodes[bl2].neighbors[1] = tl1 + 1
                nodes[tr1].neighbors[1] = br2 + 1
                nodes[br2].neighbors[0] = tr1 + 1


    for el in elements:  # Set element material properties
        el.E, el.nu, el.density, el.thickness, el.alfaC = E, nu, density, thickness, alfaC
        if stiffdef: el.stiffdef = stiffdef
    for nd in nodes:  # Assign load and supports
        for bx, by, val in load:
            if all((nd.x > bx[0], nd.x < bx[1], nd.y > by[0], nd.y < by[1])):
                nd.F_ext = val
        for bx, by, val in supports:
            if all((nd.x > bx[0], nd.x < bx[1], nd.y > by[0], nd.y < by[1])):
                nd.supports = val

    dom = type_dict[domtype](c_n=c_n, c_s=c_s, nodes=nodes, elements=elements)
    if spring_stiffness_factor:
        max_k = 0.
        for el in elements:
            if not hasattr(el, 'K'):
                el.set_matrices()
            max_k = max(max_k, np.max(el.K))
        c_n = max_k * spring_stiffness_factor
        c_s = max_k * spring_stiffness_factor
        dom.set_contacts(c_n, c_s)
    return dom

def convert_CDEM_to_FEM_domain(CDEM_dom):
    els = copy.deepcopy(CDEM_dom.elements)
    nds = copy.deepcopy(CDEM_dom.nodes)
    isclose = lambda n1, n2, tol = 1e-6: (((n2[0] - n1[0]) ** 2 + (n2[1] - n1[1]) ** 2) ** 0.5) < tol

    # Create a mapping of the CDEM nodes to FEM nodes (CDEM nodes with same coordinates condensate into one FEM node)
    mapping = {}
    new_nds = [nds[0]]
    for id1, nd1 in enumerate(nds):
        for id2, nd2 in enumerate(new_nds):
            if isclose((nd1.x, nd1.y), (nd2.x, nd2.y)):
                mapping[id1] = id2
                break
            else:
                if id2 == (len(new_nds) - 1):
                    new_nds.append(copy.deepcopy(nd1))
                    mapping[id1] = id2 + 1

    for el in els:
        el.nodes = [mapping[n - 1] + 1 for n in el.nodes]
    return DomainFEM(nodes=new_nds, elements=els)

