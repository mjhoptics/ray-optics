#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" First order paraxial design space

.. Created on Sat Mar 31 21:14:42 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np

import rayoptics.optical.model_constants as mc
import rayoptics.parax.firstorder as fo
from rayoptics.seq.gap import Gap
from rayoptics.elem.surface import Surface

from rayoptics.util.line_intersection import get_intersect
from rayoptics.util.misc_math import normalize
from numpy.linalg import norm


def bbox_from_poly(poly):
    minx, miny = np.min(poly, axis=0)
    maxx, maxy = np.max(poly, axis=0)
    return np.array([[minx, miny], [maxx, maxy]])


def calculate_slope(x, y, line_type):
    """ x=ybar, y=y  """
    if line_type == 'stop':
        k = x/y
        return k, np.array([[1, -k], [0, 1]]).T
    elif line_type == 'object_image':
        k = y/x
        return k, np.array([[1, 0], [-k, 1]]).T
    else:
        k = 0
        return 0, np.array([[1, 0], [0, 1]])


class ParaxialModel():
    def __init__(self, opt_model, opt_inv=1.0, seq_mapping=None, **kwargs):
        self.opt_model = opt_model
        self.seq_model = opt_model.seq_model
        self.seq_mapping = seq_mapping
        self.layers = {'ifcs': self} if seq_mapping is None else None
        self.sys = []
        self.ax = []
        self.pr = []
        self.opt_inv = opt_inv
        self.y_star = None
        self.ybar_star = None

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['seq_model']
        del attrs['layers']
        return attrs

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        self.seq_model = opt_model.seq_model
        if not hasattr(self, 'seq_mapping'):
            self.seq_mapping = None
        if not hasattr(self, 'layers'):
            self.layers = {'ifcs': self}
        if not hasattr(self, 'y_star'):
            self.y_star, self.ybar_star = self.calc_object_and_pupil(0)

    def update_model(self, **kwargs):
        src_model = kwargs.get('src_model', None)
        num_ifcs = len(self.seq_model.ifcs)
        if (num_ifcs > 2 and 
            (len(self.sys) != num_ifcs or src_model is not self)):
            self.build_lens()
        if kwargs.get('build', None) != 'update':
            self.layers = {'ifcs': self}

    def sync_to_seq(self, seq_model):
        self.build_lens()
        
    def set_from_specsheet(self, specsheet):
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr

        n_0 = sys[0][mc.indx]
        n_k = sys[-1][mc.indx]
        thi_0 = -self.y_star*self.ybar_star/self.opt_inv
        yu, yu_bar = specsheet.get_parax_start_data(thi_0, n_0, n_k)

        ax_ray[0][mc.ht] = y0 = yu[mc.ht]
        ax_ray[0][mc.slp] = nu0 = n_0*yu[mc.slp]
        pr_ray[0][mc.ht] = y_bar0 = yu_bar[mc.ht]
        pr_ray[0][mc.slp] = nu_bar0 = n_0*yu_bar[mc.slp]
        self.opt_inv = nu_bar0*y0 - nu0*y_bar0
        self.paraxial_trace()

        self.y_star, self.ybar_star = self.calc_object_and_pupil(0)

    def build_lens(self):
        # rebuild the `sys` description from the seq_model path
        self.sys = sys = self.seq_path_to_paraxial_lens(self.seq_model.path())

        # precalculate the reduced forms of the paraxial axial and chief rays
        parax_data = self.opt_model['analysis_results']['parax_data']
        if parax_data is not None:
            ax_ray, pr_ray, fod = parax_data
            self.opt_inv = fod.opt_inv
    
            self.ax = []
            self.pr = []
            for i in range(0, len(sys)):
                n = sys[i][mc.indx]
                self.ax.append([ax_ray[i][mc.ht], n*ax_ray[i][mc.slp]])
                self.pr.append([pr_ray[i][mc.ht], n*pr_ray[i][mc.slp]])
            self.y_star, self.ybar_star = self.calc_object_and_pupil(0)

    def init_from_nodes(self, node_list, type_sel, rndx_and_imode=None):
        """Construct a diagram using `nodes`, a list of diagram vertices. """
        node_list = np.array(node_list)
        if len(node_list.shape) == 1 or len(node_list.shape) == 3:
            nodes = node_list[type_sel]
        else:
            nodes = node_list
        max_nodes = len(nodes)
        ax = [[0., 0.] for i in range(max_nodes)]
        pr = [[0., 0.] for i in range(max_nodes)]
        sys = [[0., 0., 1., 'transmit'] for i in range(max_nodes)]
        if rndx_and_imode is None:
            rndx_and_imode = [(1., 'transmit')  for i in range(max_nodes)]
        for i, vni in enumerate(zip(nodes, rndx_and_imode)):
            vertex, ni = vni
            ax[i][type_sel] = vertex[1]
            pr[i][type_sel] = vertex[0]
            sys[i] = [0.0, 0.0, *ni]

        # adjust object distance and object node
        slp0 = ((ax[1][mc.ht]-ax[0][mc.ht])/
                (pr[1][mc.ht]-pr[0][mc.ht]))
        inc_ybar0 = ax[0][mc.ht]/slp0
        ax[0][mc.ht] = 0.0
        pr[0][mc.ht] = pr[0][mc.ht] - inc_ybar0
        nodes[0] = np.array([pr[0][mc.ht], ax[0][mc.ht]])

        self.ax = ax
        self.pr = pr
        self.sys = sys

        self.nodes_to_parax(nodes, type_sel)
        self.sys[0][mc.rmd] = 'dummy'
        self.sys[-1][mc.rmd] = 'dummy'

        return self

    def get_num_nodes(self):
        return len(self.ax)

    def parax_to_nodes(self, type_sel):
        """ render the paraxial model into a node list """
        nodes = [[x[type_sel], y[type_sel]] for x, y in zip(self.pr, self.ax)]
        nodes = np.array(nodes)
        return nodes

    def replace_node_with_dgm(self, e_node, dgm, **kwargs):
        """ Update the parax model from the node list, `nodes`. """
        prx, dgm_pkg = dgm
        pm, node_idx, type_sel = prx
        node_list, sys_data = dgm_pkg
        node_list = np.array(node_list)
        if len(node_list.shape) == 1 or len(node_list.shape) == 3:
            nodes = node_list[type_sel]
        else:
            nodes = node_list

        self.delete_node(node_idx)

        for i, ns in enumerate(zip(nodes, sys_data), start=node_idx):
            node, sys = ns
            self.add_node(i, node, type_sel, sys)

        return self

    def nodes_to_parax(self, node_list, type_sel):
        """ Update the parax model from the node list, `nodes`. """
        node_list = np.array(node_list)
        if len(node_list.shape) == 1 or len(node_list.shape) == 3:
            nodes = node_list[type_sel]
        else:
            nodes = node_list
            
        if type_sel == mc.ht:
            for i, node in enumerate(nodes):
                self.apply_ht_dgm_data(i, node)
        elif type_sel == mc.slp:
            ht_nodes = node_list[mc.ht]
            self.pr[0][mc.ht] = ht_nodes[0][0]
            self.ax[0][mc.ht] = ht_nodes[0][1]            
            for i, node in enumerate(nodes):
                self.apply_slope_dgm_data(i, node)
            self.pr[-1][mc.ht] = ht_nodes[-1][0]
            self.ax[-1][mc.ht] = ht_nodes[-1][1]
            self.sys[-2][mc.tau] = (ht_nodes[-2][1]*ht_nodes[-1][0] -
                                    ht_nodes[-1][1]*ht_nodes[-2][0])/self.opt_inv
        return self

    def get_pt(self, idx, type_sel=mc.ht):
        return [self.pr[idx][type_sel], self.ax[idx][type_sel]]

    def set_pt(self, idx, pt, type_sel=mc.ht):
        self.pr[idx][type_sel] = pt[0]
        self.ax[idx][type_sel] = pt[1]

    def apply_scale_factor(self, scale_factor):
        """ Apply scale factor to height diagram and update parax_model. """
        type_sel = mc.ht
        nodes = self.parax_to_nodes(type_sel)
        nodes = scale_factor*nodes
        self.opt_inv *= scale_factor
        self.nodes_to_parax(nodes, type_sel)

    def get_gap_for_node(self, node_idx, dgm_type):
        if self.seq_mapping is None:
            # just return the node in the seq_model
            gap_idx = node_idx
        else:  # use the node_defs to map to the seq_model
            node_defs, map_to_ifcs = self.seq_mapping
            kernel = node_defs[node_idx]
            if len(kernel) == 1:
                gap_idx = kernel[0]
            elif len(kernel) == 2:
                prev_gap_idx, after_gap_idx = kernel
                gap_idx = prev_gap_idx if dgm_type == 'ht' else after_gap_idx
            elif len(kernel) == 3:
                idx, prev_gap_idx, after_gap_idx = kernel
                gap_idx = idx if idx != after_gap_idx else after_gap_idx
        return self.seq_model.gaps[gap_idx], self.seq_model.z_dir[gap_idx]

    # --- add/delete points from diagram
    def add_node(self, node, new_vertex, type_sel, sys_data):
        """ Add a node in the paraxial data structures """

        self.sys.insert(node, [0.0, 0.0, *sys_data])

        ax_node = [0.0, 0.0]
        ax_node[type_sel] = new_vertex[1]
        self.ax.insert(node, ax_node)

        pr_node = [0.0, 0.0]
        pr_node[type_sel] = new_vertex[0]
        self.pr.insert(node, pr_node)

        if type_sel == mc.ht:
            self.apply_ht_dgm_data(node, new_vertex=new_vertex)
        elif type_sel == mc.slp:
            self.apply_slope_dgm_data(node, new_vertex=new_vertex)

        if sys_data[1] == 'reflect':
            for i in range(node, len(self.sys)):
                self.sys[i][mc.indx] = -self.sys[i][mc.indx]

        return node

    def assign_object_to_node(self, node, type_sel, factory, **inputs):
        """ create a new element from `factory` and replace `node` with it """

        # extract optical properties of node
        n = self.sys[node][mc.indx]
        power = self.sys[node][mc.pwr]
        thi = n*self.sys[node][mc.tau]
        sd = abs(self.ax[node][mc.ht]) + abs(self.pr[node][mc.ht])

        # create an element with the node's properties
        descriptor = factory(power=power, sd=sd, prx=(self, node, type_sel))
        seq, elm, e_nodez, dgm = descriptor

        # insert the path sequence and elements into the
        #  sequential and element models
        kwargs = {'idx': node, 't': thi, **inputs, 'src_model': self}
        self.opt_model.insert_ifc_gp_ele(*descriptor, **kwargs)

        self.compute_signed_rindx()
            
        n_before = self.sys[node-1][mc.indx]
        thi_before = n_before*self.sys[node-1][mc.tau]
        self.seq_model.gaps[node-1].thi = thi_before

        seq_end = node + len(seq) - 1
        n_following = self.sys[seq_end][mc.indx]
        thi_following = n_following*self.sys[seq_end][mc.tau]
        self.seq_model.gaps[seq_end].thi = thi_following

        return descriptor, kwargs

    def compute_signed_rindx(self):
        """Reset the state of the refractive index array.
        
        This method resets the signs of the refractive indices so that they are
        negative following an odd number of reflections, but positive otherwise.
        """
        flip = 1
        for ss in self.sys:
            ss[mc.indx] = abs(ss[mc.indx])
            if ss[mc.rmd] == 'reflect':
                flip = -flip
            if flip < 0:
                ss[mc.indx] = -ss[mc.indx]

    def replace_node_with_seq(self, node, sys_seq, pp_info):
        """ replaces the data at node with sys_seq """
        sys = self.sys
        ax = self.ax
        pr = self.pr
        if len(sys_seq) == 1:
            sys[node] = sys_seq[0]
        else:
            opt_inv = self.opt_inv
            efl, pp1, ppk, ffl, bfl = pp_info[2]
            sys[node-1][mc.tau] -= pp1/sys[node-1][mc.indx]
    
            # sys_seq[-1][tau] = sys[node][tau] - ppk/sys_seq[-1][indx]
            p0 = [ax[node][mc.ht], pr[node][mc.ht]]
            pn = [ax[node+1][mc.ht], pr[node+1][mc.ht]]
            slp0 = [ax[node][mc.slp], pr[node][mc.slp]]
    
            for n, ss in enumerate(sys_seq[:-1], start=node):
                sys.insert(n, ss)
    
                ax_ht = ax[n-1][mc.ht] + sys[n-1][mc.tau]*ax[n-1][mc.slp]
                ax_slp = ax[n-1][mc.slp] - ax_ht*sys[n][mc.pwr]
                ax.insert(n, [ax_ht, ax_slp])
    
                pr_ht = pr[n-1][mc.ht] + sys[n-1][mc.tau]*pr[n-1][mc.slp]
                pr_slp = pr[n-1][mc.slp] - pr_ht*sys[n][mc.pwr]
                pr.insert(n, [pr_ht, pr_slp])
    
            # replace the original node data
            ax[n+1][mc.slp] = slp0[0]
            pr[n+1][mc.slp] = slp0[1]
            sys[n+1][mc.pwr] = (ax[n][mc.slp]*pr[n+1][mc.slp] -
                             ax[n+1][mc.slp]*pr[n][mc.slp])/opt_inv
            # sys_seq[-1][pwr]
            p1 = [ax[n][mc.ht], pr[n][mc.ht]]
            p2 = [ax[n][mc.ht]+ax[n][mc.slp], pr[n][mc.ht]+pr[n][mc.slp]]
            p2int = np.array(get_intersect(p1, p2, p0, pn))
            ax[n+1][mc.ht] = p2int[0]
            pr[n+1][mc.ht] = p2int[1]
            sys[n][mc.tau] = (
                (ax[n][mc.ht]*pr[n+1][mc.ht] - ax[n+1][mc.ht]*pr[n][mc.ht])/opt_inv)
            sys[n+1][mc.tau] = (p2int[0]*pn[1] - pn[0]*p2int[1])/opt_inv

    def get_object_for_node(self, node):
        ''' basic 1:1 relationship between seq and parax model sequences '''
        ifc = self.seq_model.ifcs[node]
        e_node = self.opt_model.part_tree.parent_node(ifc)
        args = [[ifc, None, None, 1, 1]], [e_node.id], e_node
        kwargs = {'idx': node}
        return args, kwargs

    def delete_node(self, surf):
        """ delete the node at position surf """
        del self.sys[surf]
        del self.ax[surf]
        del self.pr[surf]

    # --- edit diagram points
    def apply_ht_dgm_data(self, surf, new_vertex=None):
        """ This routine calculates all data dependent on the input
            height coordinates (y,ybar) at surface surf.
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr
        opt_inv = self.opt_inv
        ht = mc.ht
        tau = mc.tau
        slp = mc.slp
        pwr = mc.pwr

        if new_vertex is not None:
            pr_ray[surf][ht] = new_vertex[0]
            ax_ray[surf][ht] = new_vertex[1]

        nsm1 = len(sys) - 1
        if surf == 0:
            surf += 1

        p = surf - 1
        c = surf

        sys[p][tau] = ((ax_ray[p][ht]*pr_ray[c][ht] -
                        ax_ray[c][ht]*pr_ray[p][ht]) / opt_inv)
        ax_ray[p][slp] = (ax_ray[c][ht] - ax_ray[p][ht])/sys[p][tau]
        pr_ray[p][slp] = (pr_ray[c][ht] - pr_ray[p][ht])/sys[p][tau]

        if (surf > 1):
            p2 = surf - 2
            sys[p][pwr] = ((ax_ray[p2][slp]*pr_ray[p][slp] -
                            ax_ray[p][slp]*pr_ray[p2][slp])
                           / opt_inv)

        if (surf < nsm1):
            s = surf + 1

            sys[c][tau] = (ax_ray[c][ht]*pr_ray[s][ht] -
                           ax_ray[s][ht]*pr_ray[c][ht])/opt_inv
            ax_ray[c][slp] = (ax_ray[s][ht] - ax_ray[c][ht])/sys[c][tau]
            pr_ray[c][slp] = (pr_ray[s][ht] - pr_ray[c][ht])/sys[c][tau]
            sys[c][pwr] = (ax_ray[p][slp]*pr_ray[c][slp] -
                           ax_ray[c][slp]*pr_ray[p][slp])/opt_inv

            if s < nsm1:
                sys[s][pwr] = (ax_ray[c][slp]*pr_ray[s][slp] -
                               ax_ray[s][slp]*pr_ray[c][slp])/opt_inv

        else:
            ax_ray[c][slp] = ax_ray[p][slp]
            pr_ray[c][slp] = pr_ray[p][slp]
            sys[c][pwr] = 0
            sys[c][tau] = 0

    def apply_slope_dgm_data(self, surf, new_vertex=None):
        """ This routine calculates all data dependent on the input
            slope coordinates (nu,nubar) at surface surf.
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr
        opt_inv = self.opt_inv
        ht = mc.ht
        tau = mc.tau
        slp = mc.slp
        pwr = mc.pwr

        if new_vertex is not None:
            pr_ray[surf][slp] = new_vertex[0]
            ax_ray[surf][slp] = new_vertex[1]

        nsm1 = len(sys) - 1
        if nsm1 == 0:
            p = 0
            c = 1
            ax_ray[c][ht] = ax_ray[p][slp]*sys[p][tau] + ax_ray[p][ht]
            pr_ray[c][ht] = pr_ray[p][slp]*sys[p][tau] + pr_ray[p][ht]

        elif (surf < nsm1):
            if (surf == 0):
                surf += 1

            p = surf - 1
            c = surf

            sys[c][pwr] = (ax_ray[p][slp]*pr_ray[c][slp] -
                           ax_ray[c][slp]*pr_ray[p][slp])/opt_inv
            ax_ray[c][ht] = (ax_ray[p][slp] - ax_ray[c][slp])/sys[c][pwr]
            pr_ray[c][ht] = (pr_ray[p][slp] - pr_ray[c][slp])/sys[c][pwr]
            sys[p][tau] = (ax_ray[p][ht]*pr_ray[c][ht] -
                           ax_ray[c][ht]*pr_ray[p][ht])/opt_inv
           
            s = surf + 1
            ax_ray[s][ht] = ax_ray[c][slp]*sys[c][tau] + ax_ray[c][ht]
            pr_ray[s][ht] = pr_ray[c][slp]*sys[c][tau] + pr_ray[c][ht]

    def process_seq_mapping(self, nodes):
        """The composite `nodes` are mapped and applied to the ifcs layer. """
        node_defs, map_to_ifcs = self.seq_mapping
        nodes_ifcs = []
        for i, kernel in enumerate(map_to_ifcs):
            idx, nidx, t1, t2 = kernel
            if t1 == 0:
                new_node = t2*np.array(nodes[nidx])
            else:
                l1_pt1 = np.array(nodes[nidx])
                l1_pt2 = np.array(nodes[nidx+1])
                new_intersection_pt = t1*(l1_pt2 - l1_pt1) + l1_pt1
                new_node = t2*new_intersection_pt
            nodes_ifcs.append(new_node)

        # print(f'nodes_ifcs {nodes_ifcs}')
        seq_nodes = [nodes_ifcs, []]
        self.opt_model['parax_model'].nodes_to_parax(seq_nodes, mc.ht)

    def update_composite_node_fct(self, type_sel):
        """ returns the appropriate 'apply' fct for the `type_sel`. """            
        if type_sel == mc.ht:
            def apply_ht_data(node, new_vertex):
                self.apply_ht_dgm_data(node, new_vertex)
                if self.seq_mapping is not None:
                    ht_nodes = self.parax_to_nodes(mc.ht)
                    self.process_seq_mapping(ht_nodes)
            return apply_ht_data

        elif type_sel == mc.slp:
            def apply_slp_data(node, new_vertex):
                # editing only allowed for interface layer currently
                self.apply_slope_dgm_data(node, new_vertex)
                if self.seq_mapping is not None:
                    ht_nodes = self.parax_to_nodes(mc.ht)
                    self.process_seq_mapping(ht_nodes)
            return apply_slp_data

    def update_rindex(self, surf):
        """Update the refractive index using the `gap` at *surf*."""
        gap = self.seq_model.gaps[surf]
        wvl = self.seq_model.central_wavelength()
        self.sys[surf][mc.indx] = gap.medium.rindex(wvl)

    # ParaxTrace() - This routine performs a paraxial raytrace from object
    #                (surface 0) to image.  The last operation is a
    #                transfer to the image surface.
    def paraxial_trace(self):
        """ regenerate paraxial axial and chief rays from power and reduced
            distance
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr
        ht = mc.ht
        tau = mc.tau
        slp = mc.slp
        pwr = mc.pwr

        nsm1 = len(sys) - 1

        # Transfer from object
        c = 0
        s = 1
        ax_ray[s][ht] = ax_ray[c][ht] + sys[c][tau]*ax_ray[c][slp]
        pr_ray[s][ht] = pr_ray[c][ht] + sys[c][tau]*pr_ray[c][slp]

        for i in range(1, len(sys)):

            p = c
            c = s

            # Refraction
            ax_ray[c][slp] = ax_ray[p][slp] - ax_ray[c][ht]*sys[c][pwr]
            pr_ray[c][slp] = pr_ray[p][slp] - pr_ray[c][ht]*sys[c][pwr]

            # Transfer
            if (i < nsm1):
                s += 1
                ax_ray[s][ht] = ax_ray[c][ht] + sys[c][tau]*ax_ray[c][slp]
                pr_ray[s][ht] = pr_ray[c][ht] + sys[c][tau]*pr_ray[c][slp]

    # --- list output
    def list_model(self):
        """ list the paraxial axial and chief rays, and power, reduced distance
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr

        header_str = (
            "           ax_ht        pr_ht       ax_slp       pr_slp"
            "         power          tau        index    type")

        print(header_str)
        for i, aps in enumerate(zip(ax_ray, pr_ray, sys)):
            ax, pr, sy = aps
            print(f"{i:2}: {ax[mc.ht]:12.6g} {pr[mc.ht]:12.6g} "
                  f"{ax[mc.slp]:12.6g} {pr[mc.slp]:12.6g} "
                  f"{sy[mc.pwr]:13.7g} {sy[mc.tau]:12.6g} "
                  f"{sy[mc.indx]:12.5f}    {sy[mc.rmd]}")

    def list_lens(self):
        """ list the paraxial axial and chief rays, and power, reduced distance
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr
        ht = mc.ht
        tau = mc.tau
        slp = mc.slp
        pwr = mc.pwr
        indx = mc.indx
        rmd = mc.rmd

        print("       ax_ray_ht    ax_ray_slp")
        for i in range(0, len(ax_ray)):
            print("{:2}: {:12.6g}  {:12.6g}".format(i, ax_ray[i][ht],
                                                    ax_ray[i][slp]))
        
        print("\n       pr_ray_ht    pr_ray_slp")
        for i in range(0, len(pr_ray)):
            print("{:2}: {:12.6g}  {:12.6g}".format(i, pr_ray[i][ht],
                                                    pr_ray[i][slp]))

        print("\n            power           tau        index    type")
        for i in range(0, len(sys)):
            print("{:2}: {:13.7g}  {:12.5g} {:12.5f}    {}".format(
                i, sys[i][pwr], sys[i][tau], sys[i][indx], sys[i][rmd]))

    def list_sys_seq(self):
        print("             c            t        medium     mode")
        fmt = "{:2}: {:12.7g} {:#12.6g} {:12.5f} {:>10s}"
        sys = self.sys
        n_before = sys[0][mc.indx]
        for i, sg in enumerate(self.sys):
            n_after = sys[i][mc.indx]
            thi = n_after*sys[i][mc.tau]
            delta_n = n_after - n_before
            cv = sys[i][mc.pwr] if delta_n == 0 else sys[i][mc.pwr]/delta_n
            print(fmt.format(i, cv, thi, n_after, sys[i][mc.rmd]))
            n_before = n_after

    def first_order_data(self):
        """List out the first order imaging properties of the model."""
        self.opt_model['analysis_results']['parax_data'].fod.list_first_order_data()

    # --- convert to/from sequential model

    def seq_path_to_paraxial_lens(self, path):
        """ returns lists of power, reduced thickness, signed index and refract
            mode
        """
        sys = []
        for i, sg in enumerate(path):
            ifc, gap, _, rndx, z_dir = sg
            imode = ifc.interact_mode
            power = ifc.optical_power
            if gap:
                n_after = rndx if z_dir > 0 else -rndx
                tau = gap.thi/n_after
                sys.append([power, tau, n_after, imode])
            else:
                sys.append([power, 0.0, sys[-1][mc.indx], imode])
        return sys

    def paraxial_lens_to_seq_model(self):
        """ Applies a paraxial lens spec (power, reduced distance) to the model
        """
        sys = self.sys
        ax_ray = self.ax
        n_before = sys[0][mc.indx]
        slp_before = self.ax[0][mc.slp]
        for i, sg in enumerate(self.seq_model.path()):
            ifc, gap, _, _, _ = sg
            if gap:
                n_after = sys[i][mc.indx]
                slp_after = ax_ray[i][mc.slp]
                gap.thi = n_after*sys[i][mc.tau]

                ifc.set_optical_power(sys[i][mc.pwr], n_before, n_after)
                ifc.from_first_order(slp_before, slp_after, ax_ray[i][mc.ht])

                n_before = n_after
                slp_before = slp_after

    def lens_from_dgm(self, node, bending=0., th=None, sd=1., **kwargs):
        """ Return single lens constructional parameters from parax_model. 
        
        This method uses a method for thickening a lens element developed by
        Lopez-Lopez in his PhD dissertation on the |yybar| diagram. 
        
        The variables are named following the notation in the thesis. He used
        `z` for the |yybar| coordinates; the reduced slopes and system
        parameters include a (inverse) factor of the optical invariant, and
        are capitalized here: T, W, Pwr
        """
        from opticalglass.modelglass import ModelGlass  # type: ignore
        if 'med' in kwargs:
            med = kwargs['med']
            if med is None:
                med = ModelGlass(1.517, 64.2, '517642')
            else:
                from rayoptics.seq.medium import decode_medium
                med = decode_medium(med)
        else:
            med = ModelGlass(1.517, 64.2, '517642')
        
        nom_wvl = kwargs.get('nom_wvl', 'd')
        rndx = med.rindex(nom_wvl)

        opt_inv = self.opt_inv
        z0 = np.array(self.get_pt(node-1, type_sel=mc.ht))
        z1 = np.array(self.get_pt(node, type_sel=mc.ht))
        z2 = np.array(self.get_pt(node+1, type_sel=mc.ht))
    
        To = np.cross(z1, z0)
        Wo = (z1 - z0)/To

        Ti = np.cross(z2, z1)
        Wi = (z2 - z1)/Ti
        Pwr = np.cross(Wi, Wo)
        # print(f"To={To:8.4f}, Ti={Ti:8.4f}, opt_inv={opt_inv:8.4f}")
        # print(f"Pwr={Pwr:8.6f}, efl={1/(Pwr*opt_inv):8.4f}")

        zf = -Wi/Pwr
        zfp = Wo/Pwr
        # print(f"zf={zf}, zfp={zfp}")

        sd = kwargs.get('sd', norm(z1))
        th = kwargs.get('th', None)
        if th is None:
            th = sd/5

        X = bending
        T_lens = (th/rndx)*opt_inv
        p = T_lens*Pwr
        Q = 1 + p*(X**2 - 1)
        # brkt = 1/(1 + np.sqrt(Q))
        # print(f"X={X:8.4f}, T_lens={T_lens:8.4f}, p={p:8.4f}, Q={Q:8.4f}, brkt={brkt:8.4f}")

        FWo = (X + 1)/(1 + np.sqrt(Q))
        FWi = -(X - 1)/(1 + np.sqrt(Q))
        FZo = (X + np.sqrt(Q))/(X + 1)
        FZi = (X - np.sqrt(Q))/(X - 1)
        # print(f"FWo={FWo:8.4f}, FWi={FWi:8.4f}, "
        #       f"FZo={FZo:8.4f}, FZi={FZi:8.4f}")

        Wx = FWo*Wi + FWi*Wo
        z10 = zf + FZo*zfp
        z11 = FZi*zf + zfp
        P10 = FWo*Pwr
        P11 = FWi*Pwr
        # To10 = np.cross(z10, z0)
        # T11i = np.cross(z2, z11)
        
        # pp1 = (To - To10)/opt_inv
        # ppk = (Ti - T11i)/opt_inv
        # print(f"z0={z0}, z10={z10}, z1={z1}, z11={z11}, z2={z2}")
        # print(f"To10={To10:8.4f}, T11i={T11i:8.4f}, "
        #       f"pp1={pp1:8.4f}, ppk={ppk:8.4f}")
    
        cv1 = P10*opt_inv/(rndx-1)
        cv2 = P11*opt_inv/(1-rndx)
        thi = rndx*T_lens/opt_inv
        lens = cv1, cv2, thi, rndx, sd
        # print(f"cv1={cv1:8.6f}, cv2={cv2:8.6f}, thi={thi:8.4f}, "
        #       f"rndx={rndx:6.4f}, sd={sd:8.4f}")
        
        dgm_pkg = [z10, z11], [[rndx, 'transmit'], [1.0, 'transmit']]
        return dgm_pkg, lens

    def calc_object_and_pupil(self, idx):
        """ calculate axial intercepts of line between idx and idx+1 """
        k = ((self.ax[idx+1][mc.ht]-self.ax[idx][mc.ht])/
             (self.pr[idx+1][mc.ht]-self.pr[idx][mc.ht]))
        y_star = self.ax[idx][mc.ht] - k*self.pr[idx][mc.ht]
        ybar_star = self.pr[idx][mc.ht] - self.ax[idx][mc.ht]/k
        # print(f"{idx}: y_star={y_star}, ybar_star={ybar_star}, k={k}")
        return y_star, ybar_star

    # --- power and thickness solves
    def pwr_slope_solve(self, ray, surf, slp_new):
        p = ray[surf-1]
        c = ray[surf]
        pwr = (p[mc.slp] - slp_new)/c[mc.ht]
        return pwr

    def pwr_ht_solve(self, ray, surf, ht_new):
        sys = self.sys
        p = ray[surf-1]
        c = ray[surf]
        slp_new = (ht_new - c[mc.ht])/sys[surf][mc.tau]
        pwr = (p[mc.slp] - slp_new)/ht_new
        return pwr

    def thi_ht_solve(self, ray, surf, ht_new):
        c = ray[surf]
        thi = (ht_new - c[mc.ht])/c[mc.slp]
        return thi

    # --- calculations
    def compute_principle_points(self, seq):
        """ Returns paraxial p and q rays, plus partial first order data.

        Args:
            seq: a sequence containing interfaces and gaps to be traced.
                  for each iteration, the sequence should return a
                  list containing: **Intfc, Gap, Trfm, Index, Z_Dir**

        Returns:
            (p_ray, q_ray, (efl, pp1, ppk, ffl, bfl))

            - p_ray: [ht, slp, aoi], [1, 0, -]
            - q_ray: [ht, slp, aoi], [0, 1, -]
            - efl: effective focal length
            - pp1: distance of front principle plane from 1st interface
            - ppk: distance of rear principle plane from last interface
            - ffl: front focal length
            - bfl: back focal length
        """
        n_0 = seq[0][mc.Indx]
        z_dir_before = seq[0][mc.Zdir]
        n_k = seq[-1][mc.Indx]
        z_dir_k = seq[-1][mc.Zdir]
        path = [[Surface(), Gap(), None, n_0, z_dir_before]]
        path.extend(seq[1:])
        pp_info = fo.compute_principle_points(iter(path),
                                              n_0=z_dir_before*n_0,
                                              n_k=z_dir_k*n_k)
        pp = pp_info[2]
        print(f"efl={pp[0]:8.4f}, ffl={pp[3]:8.4f}, bfl={pp[4]:8.4f}, "
              f"pp1={pp[1]:8.4f}, ppk={pp[2]:8.4f}")
        return pp_info

    def compute_principle_points_from_dgm(self, os_idx=0, is_idx=-2):
        """ Calculate focal points and principal points.

        Args:
            os_idx: object space gap index
            is_idx: image space gap index

        Returns:
            (efl, pp1, ppk, ffl, bfl)

            - efl: effective focal length
            - pp1: distance of front principle plane from 1st interface
            - ppk: distance of rear principle plane from last interface
            - ffl: front focal length
            - bfl: back focal length
        """
        n_o = self.sys[os_idx][mc.indx]
        n_i = self.sys[is_idx][mc.indx]

        opt_inv = self.opt_inv

        z0 = np.array(self.get_pt(os_idx))
        z1o = np.array(self.get_pt(os_idx+1))
        z1i = np.array(self.get_pt(is_idx))
        z2 = np.array(self.get_pt(is_idx+1))

        z1 = np.array(get_intersect(z0, z1o, z1i, z2))

        To = np.cross(z1, z1o)
        Wo = (z1 - z1o)/To

        Ti = np.cross(z1i, z1)
        Wi = (z1i - z1)/Ti
        Pwr = np.cross(Wi, Wo)
        
        zfo = -Wi/Pwr
        zfi = Wo/Pwr

        power = Pwr*opt_inv
        pp1 = n_o*To/opt_inv
        ppk = n_i*Ti/opt_inv
        efl = 1/power
        ffl = pp1 - efl
        bfl = efl - ppk
        pp_info = efl, pp1, ppk, ffl, bfl
        # print(f"efl={pp_info[0]:8.4f}, ffl={pp_info[3]:8.4f}, bfl={pp_info[4]:8.4f}, "
        #       f"pp1={pp_info[1]:8.4f}, ppk={pp_info[2]:8.4f}")
        return pp_info


    def apply_conjugate_shift(self, nodes, k, mat, line_type):
        sheared_nodes = self.calc_conjugate_shift(nodes, k, mat, line_type)

        if self.seq_mapping is not None:
            self.process_seq_mapping(sheared_nodes)

        parax_model = self.opt_model['parax_model']
        parax_model.paraxial_lens_to_seq_model()
        self.opt_model['optical_spec'].sync_to_parax(parax_model)

    def calc_conjugate_shift(self, nodes, k, mat, line_type):
        sys = self.sys
        ax = self.ax
        pr = self.pr
        opt_inv = self.opt_inv

        # apply the shear transformation to the original shape
        sheared_nodes = np.matmul(nodes, mat)

        # for an object shift, object and image distances change,
        #  calculate them before revising the ray slope values
        if line_type == 'object_image':
            conj_line = np.array([0., 0.]), np.array([1., k])

            # update object distance and object y and ybar
            sheared_obj = get_intersect(nodes[0], nodes[1], *conj_line)
            pr[0][mc.ht], ax[0][mc.ht] = sheared_obj[0], 0.
            sys[0][mc.tau] = ((-sheared_nodes[1][1]*sheared_obj[0]) / opt_inv)

            # update image distance and image y and ybar
            sheared_img = get_intersect(nodes[-2], nodes[-1], *conj_line)
            pr[-1][mc.ht], ax[-1][mc.ht] = sheared_img[0], 0.
            sys[-2][mc.tau] = ((sheared_nodes[-2][1]*sheared_img[0]) / opt_inv)

        else:  # pupil shift, update object values here
            pr[0][mc.ht], ax[0][mc.ht] = sheared_nodes[0]

        # update height and slope values over the diagram
        for i, sheared_node in enumerate(sheared_nodes[1:-1], start=1):
            pr[i][mc.ht], ax[i][mc.ht] = sheared_node
            pr[i-1][mc.slp] = (
                (pr[i][mc.ht] - pr[i-1][mc.ht]) / sys[i-1][mc.tau])
            ax[i-1][mc.slp] = (
                (ax[i][mc.ht] - ax[i-1][mc.ht]) / sys[i-1][mc.tau])

        # update final slope values
        pr[i][mc.slp] = (pr[i+1][mc.ht] - pr[i][mc.ht]) / sys[i][mc.tau]
        ax[i][mc.slp] = (ax[i+1][mc.ht] - ax[i][mc.ht]) / sys[i][mc.tau]
        pr[-1][mc.slp] = pr[-2][mc.slp]
        ax[-1][mc.slp] = ax[-2][mc.slp]
        
        return sheared_nodes

    def match_pupil_and_conj(self, prx, **kwargs):
        pm, node, type_sel = prx
        opt_inv = pm.opt_inv
        z0 = np.array(pm.get_pt(node-1, type_sel=mc.ht))
        z1 = np.array(pm.get_pt(node, type_sel=mc.ht))
        z2 = np.array(pm.get_pt(node+1, type_sel=mc.ht))
        To = np.cross(z1, z0)
        Wo = (z1 - z0)/To

        pp_info = self.compute_principle_points_from_dgm()
        efl, pp1, ppk, ffl, bfl = pp_info
        T_pp1 = pp1 * opt_inv
        To_new = To - T_pp1
        z1_new = z0 + To_new * Wo

        T_ppk = ppk * opt_inv
        To_new = To - T_pp1
        z1_new = z0 + To_new * Wo

        nodes = self.parax_to_nodes(type_sel)
        shear = nodes[1] - z1_new

        k_oi, mat_oi = calculate_slope(nodes[1][0], shear[1], 'object_image')
        obj_shft = self.calc_conjugate_shift(nodes, k_oi, mat_oi, 
                                             'object_image')

        k_s, mat_s = calculate_slope(shear[0], obj_shft[1][1], 'stop')
        node_list = self.calc_conjugate_shift(obj_shft, k_s, mat_s, 'stop')

        sys_data = [[sys_i[mc.indx], sys_i[mc.indx]] for sys_i in self.sys]

        dgm_pkg = node_list[1:-1], sys_data[1:-1]
        dgm = prx, dgm_pkg
        return dgm

    def paraxial_vignetting(self, rel_fov=1):
        """Calculate the vignetting factors using paraxial optics. """
        sm = self.seq_model
        min_vly = 1, None
        min_vuy = 1, None
        for i, ifc in enumerate(sm.ifcs[:-1]):
            if self.ax[i][mc.ht] != 0:
                max_ap = ifc.surface_od()
                y = self.ax[i][mc.ht]
                ybar = rel_fov * self.pr[i][mc.ht]
                ratio = (max_ap - abs(ybar))/abs(y)
                if ratio > 0:
                    if ybar < 0:
                        if ratio < min_vly[0]:
                            min_vly = ratio, i
                    elif ybar > 0:
                        if ratio < min_vuy[0]:
                            min_vuy = ratio, i
                    else:  # ybar == 0
                        if ratio < min_vly[0]:
                            min_vly = ratio, i
                        if ratio < min_vuy[0]:
                            min_vuy = ratio, i
                        
                # print(f'{i:2d}: {ratio:8.3f} {ybar:8.3f} {y:8.3f} {max_ap:8.3f}')

        # print("min_vly:", min_vly, "min_vuy:", min_vuy)
        return min_vly, min_vuy


def nodes_to_new_model(*args, **inputs):
    """ Sketch a 2d polyline and populate an opt_model from it. """
    opm = inputs['opt_model']
    pm = opm['parax_model']
    fig = inputs['figure']
    factory = inputs['factory']
    dgm = fig.diagram
    def on_finished(poly_data, line_artist):
        pm.init_from_nodes(poly_data, dgm.type_sel)
        create_seq_from_parax(opm, factory, dgm.type_sel)
        opm['optical_spec'].sync_to_parax(pm)
        fig.refresh_gui(build='rebuild', src_model=pm)
    gui_fct = inputs['gui_fct']
    gui_fct(figure=fig, do_on_finished=on_finished)


def create_seq_from_parax(opt_model, factory, type_sel):
    pm = opt_model['parax_model']
    opt_model['seq_model'].reset()
    opt_model['ele_model'].reset()
    opt_model['part_tree'].root_node.children = []
    for i in range(1, len(pm.ax)-1):
        # extract optical properties of node
        n = pm.sys[i][mc.indx]
        power = pm.sys[i][mc.pwr]
        thi = n*pm.sys[i][mc.tau]
        sd = abs(pm.ax[i][mc.ht]) + abs(pm.pr[i][mc.ht])

        # create an element with the node's properties
        
        descriptor = factory(power=power, sd=sd, prx=(pm, i, type_sel))
        seq, ele, e_node, dgm = descriptor

        n_before = pm.sys[i-1][mc.indx]
        thi_before = n_before*pm.sys[i-1][mc.tau]

        # insert the path sequence and elements into the
        #  sequential and element models
        kwargs = {'idx': i-1, 't': thi, 
                  'insert': True, 'src_model': pm, 'do_update': False}
        opt_model.insert_ifc_gp_ele(*descriptor, **kwargs)

        opt_model['seq_model'].gaps[i-1].thi = thi_before


def create_diagram_for_key(opm, key, type_sel):
    seq_mapping, nodes = generate_mapping_for_key(opm, key, type_sel)
    prx_model = build_from_yybar(opm, nodes, type_sel, seq_mapping)
    return key, prx_model    


def update_diagram_for_key(opm, key, type_sel):
    seq_mapping, nodes = generate_mapping_for_key(opm, key, type_sel)
    if key in opm['parax_model'].layers:
        prx_model = opm['parax_model'].layers[key]
        prx_model.seq_mapping = seq_mapping
        rndx_and_imode = None
        if key == 'ifcs':
            rndx_and_imode = opm['seq_model'].get_rndx_and_imode()
        prx_model.init_from_nodes(nodes, type_sel,
                                  rndx_and_imode=rndx_and_imode)
    else:
        prx_model = build_from_yybar(opm, nodes, type_sel, seq_mapping)
        opm['parax_model'].layers[key] = prx_model
    return key, prx_model    


def generate_mapping_for_key(opm, key, type_sel):
    """ generate all of the mappings, ht and slp, for *key*. """
    
    # generate the node_defs for key
    pt = opm['part_tree']
    if key == 'eles':
        air_nodes = pt.nodes_with_tag('#airgap')
        air_gap_list = [n.id.reference_idx() for n in air_nodes]
        node_defs = air_gaps_to_node_defs(air_gap_list)
    elif key == 'asm':
        air_gap_list = [n.id.reference_idx() for n in pt.root_node.children
                        if '#airgap' in n.tag]
        node_defs = air_gaps_to_node_defs(air_gap_list)
    elif key == 'sys':
        num_nodes = len(opm['seq_model'].ifcs)
        node_defs = [(0,), (0, num_nodes-2), (num_nodes-1,)]
    elif key == 'ifcs':
        nodes = [opm['parax_model'].parax_to_nodes(mc.ht), 
                 opm['parax_model'].parax_to_nodes(mc.slp)]
        return None, nodes

    # print(f'{key}:\nnode_defs_orig {node_defs}')
    node_defs, ht_nodes = get_valid_ht_nodes(opm['parax_model'], node_defs)
    # print(f'node_defs {node_defs}\nnodes {nodes}')

    map_to_ifcs = gen_ifcs_node_mapping(opm['parax_model'], node_defs, ht_nodes)
    # print(f'map_to_ifcs {map_to_ifcs}\n')

    slp_nodes = slp_nodes_from_node_defs(opm['parax_model'], node_defs)

    seq_mapping = (node_defs, map_to_ifcs)
    nodes = [np.array(ht_nodes), np.array(slp_nodes)]
    return seq_mapping, nodes


def air_gaps_to_node_defs(air_gap_list):
    """ generate the node defs for composite layers, based on airgaps. """
    prev_gap_idx = 0
    node_defs = [(prev_gap_idx,)]
    for gap_idx in air_gap_list[1:]:
        if gap_idx - prev_gap_idx < 2:
            node_defs.append((gap_idx,))
        else:
            node_defs.append((prev_gap_idx, gap_idx))
        prev_gap_idx = gap_idx
    node_defs.append((gap_idx+1,))
    return node_defs


def get_valid_ht_nodes(parax_model, node_defs_in):
    """ given the input node defs, replace non-physical thin lenses as needed."""
    node_defs = None
    node_defs_new = node_defs_in
    nodes_new = ht_nodes_from_node_defs(parax_model, node_defs_new)

    while node_defs != node_defs_new:
        node_defs = node_defs_new
        nodes = nodes_new
        node_defs_new = scan_nodes(parax_model, node_defs, nodes)
        nodes_new = ht_nodes_from_node_defs(parax_model, node_defs_new)
    return node_defs_new, nodes_new


def ht_nodes_from_node_defs(parax_model, node_defs):
    """ produce a list of ht nodes given the parax_model and node_defs.

    `node_defs` is a list of tuples, each with either one or two indices.
    if there is a single index, it is to a node in `parax_model`.
    if there are 2 indices, the first is to the gap preceding the element;
    the second is to the gap following the element (also the last interface
    of the element). The node is calculated from the intersection of the
    diagram edges corresponding to these gaps.
    
    There is no guarentee that the nodes calculated here represent a physically
    realizable system, i.e. there may be virtual airspaces.
    """
    nodes = []
    for i, kernel in enumerate(node_defs):
        if len(kernel) == 1:
            nodes.append(parax_model.get_pt(kernel[0], type_sel=mc.ht))
        elif len(kernel) == 2:
            prev_gap_idx, after_gap_idx = kernel
            l1_pt1 = parax_model.get_pt(prev_gap_idx)
            l1_pt2 = parax_model.get_pt(prev_gap_idx+1)
            l2_pt1 = parax_model.get_pt(after_gap_idx)
            l2_pt2 = parax_model.get_pt(after_gap_idx+1)
            new_node = get_intersect(l1_pt1, l1_pt2, l2_pt1, l2_pt2)
            nodes.append(new_node)
        elif len(kernel) == 3:
            idx, prev_gap_idx, after_gap_idx = kernel
            nodes.append(parax_model.get_pt(idx))
    return nodes


def slp_nodes_from_node_defs(parax_model, node_defs):
    """ produce a list of slopes given the parax_model and node_defs.

    `node_defs` is a list of tuples, each with either one or two indices.
    if there is a single index, it is to a node in `parax_model`.
    if there are 2 indices, the first is to the gap preceding the element;
    the second is to the gap following the element (also the last interface
    of the element). The node is calculated from the intersection of the
    diagram edges corresponding to these gaps.
    
    There is no guarentee that the nodes calculated here represent a physically
    realizable system, i.e. there may be virtual airspaces.
    """
    nodes = []
    for kernel in node_defs[:-1]:
        if len(kernel) == 1:
            after_gap_idx = kernel[0]
        elif len(kernel) == 2:
            prev_gap_idx, after_gap_idx = kernel
        elif len(kernel) == 3:
            idx, prev_gap_idx, after_gap_idx = kernel

        nodes.append(parax_model.get_pt(after_gap_idx, type_sel=mc.slp))

    return nodes


def scan_nodes(parax_model, node_defs, nodes):
    """scan node defs for any invalid thin elements

    Replace the first invalid thin element found with two 3 element tuples,
    signifying a thick node. The first tuple element is the index to the node
    in the `parax_model` and the last two elements are the range of indices
    in the `parax_model` covered by the thick element.
    
    Return the updated node_def list.
    """
    new_node_defs = node_defs.copy()
    xprods = np.cross(nodes[:-1], nodes[1:])
    for i, kn in enumerate(zip(node_defs, nodes)):
        kernel, node = kn
        if len(kernel) == 2:
            prev_gap_idx, after_gap_idx = kernel
            prev = 1 if prev_gap_idx == 0 else prev_gap_idx
            l1_pt1 = parax_model.get_pt(prev)
            # l1_pt2 = parax_model.get_pt(prev_gap_idx+1)
            l2_pt1 = parax_model.get_pt(after_gap_idx)
            xprod1 = np.cross(l1_pt1, l2_pt1)
            # xprodt = np.cross(l1_pt2, l2_pt1)
            # print(f'{i}: {prev}-{after_gap_idx}: ifc: {xprod1}, t: {xprodt}, thin: {xprods[i-1]}')
            if xprods[i-1] > 0 or (prev_gap_idx != 0 and
                                   2*xprod1 > xprods[i-1]):
                new_node_defs[i] = (prev_gap_idx+1,
                                    prev_gap_idx, after_gap_idx)
                new_node_defs.insert(i+1, (after_gap_idx,
                                           prev_gap_idx, after_gap_idx))
                break
    return new_node_defs


def build_from_yybar(opm, nodes, type_sel, seq_mapping):
    prx_model = ParaxialModel(opm, opt_inv=opm['pm'].opt_inv,
                              seq_mapping=seq_mapping)
    rndx_and_imode = None
    if seq_mapping is None:
        rndx_and_imode = opm['seq_model'].get_rndx_and_imode()
    prx_model.init_from_nodes(nodes, type_sel, 
                              rndx_and_imode=rndx_and_imode)
    return prx_model


def gen_ifcs_node_mapping(parax_model, node_defs, nodes):
    """Create mapping between composite diagram and interface based diagram. 
    
    Each node in the composite diagram is associated with one or a range of
    nodes in parax_model.layer['ifcs']. `node_defs` and `nodes` define
    the composite diagram. 
    
    `node_defs` is a list of tuples, one per composite node, of length 1, 2,
    or 3. The number of entries is as follows:
        
    1) the composite node maps directly to node idx in the 'ifcs' layer
    2) the composite node is generated from the previous and following gap indices
    3) the composite node is part of a thick node
    
    A **thick** node is what is used when reducing a range of interfaces to a
    single node requires virtual propagation distances. In this case, the first
    and last nodes in the range are retained in the composite diagram; interior
    nodes are scaled according to how the thick edge is stretched.
    
    Changes in the composite diagram are propagated to the underlying 'ifcs'
    layer by applying a 2D stretch to the nodes in the 'ifcs' layer. The 'ifcs'
    node is parameterized by the calculating the intersection of the composite
    edge with the vector from the origin through the composite node. 
    The scale factors are:
        
    * t1: parametric distance of the intersection point along the composite edge
    * t2: fractional distance of the composite node to the intersection point
        
    The map_to_ifcs list connects the edge in the composite diagram to the 
    'ifcs' node and the scale factors needed to update the 'ifcs' node position
    when the composite diagram changes.
    """
    def calc_scale_factors(pt1, pt2, pt_k):
        new_node = np.array(get_intersect(pt1, pt2, np.array([0., 0.]), pt_k))
        t1 = norm(new_node - pt1)/norm(pt2 - pt1)
        t2 = norm(pt_k)/norm(new_node)
        return t1, t2

    map_to_ifcs = []
    for i, kernel in enumerate(node_defs):
        if len(kernel) == 1:  # single node or mirror
            idx = kernel[0]
            map_to_ifcs.append((idx, i, 0., 1.))

        elif len(kernel) == 2:  # thin element group
            # get the defining vertices from the composite diagram
            l1_pt1 = np.array(nodes[i-1])
            l1_pt2 = np.array(nodes[i])
            l2_pt1 = np.array(nodes[i])
            l2_pt2 = np.array(nodes[i+1])

            prev_gap_idx, after_gap_idx = kernel
            for k in range(prev_gap_idx+1, after_gap_idx+1):
                pt_k = np.array(parax_model.get_pt(k, type_sel=mc.ht))
                xprod = np.cross(nodes[i], pt_k)
                if xprod >= 0.0:
                    t1, t2 = calc_scale_factors(l1_pt1, l1_pt2, pt_k)
                    map_to_ifcs.append((k, i-1, t1, t2))
                else:
                    t1, t2 = calc_scale_factors(l2_pt1, l2_pt2, pt_k)
                    map_to_ifcs.append((k, i, t1, t2))

        elif len(kernel) == 3:  # thick element group
            # A thick element group has 2 adjacent nodes. The first and last nodes
            # of the thick element correspond to the first and last nodes in the
            # interface diagram. Interface nodes between these are scaled along
            # the edge separating the 2 thick nodes. 
            map_to_ifcs.append((idx, i, 0., 1.))

            idx, prev_gap_idx, after_gap_idx = kernel
            if idx != after_gap_idx:
                # get the defining vertices from the composite diagram
                thick_pt1 = np.array(nodes[i])
                thick_pt2 = np.array(nodes[i+1])
                for k in range(idx+1, after_gap_idx):
                    pt_k = np.array(parax_model.get_pt(k, type_sel=mc.ht))
                    t1, t2 = calc_scale_factors(thick_pt1, thick_pt2, pt_k)
                    map_to_ifcs.append((k, i, t1, t2))

    return map_to_ifcs
