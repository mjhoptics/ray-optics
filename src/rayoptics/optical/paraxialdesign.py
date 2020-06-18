#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" First order paraxial design space

.. Created on Sat Mar 31 21:14:42 2018

.. codeauthor: Michael J. Hayford
"""

import numpy as np

from rayoptics.optical.model_constants import ht, slp, aoi
from rayoptics.optical.model_constants import pwr, tau, indx, rmd
import rayoptics.optical.model_constants as mc
import rayoptics.optical.firstorder as fo
from rayoptics.optical.gap import Gap
from rayoptics.optical.surface import Surface


def bbox_from_poly(poly):
    minx, miny = np.min(poly, axis=0)
    maxx, maxy = np.max(poly, axis=0)
    return np.array([[minx, miny], [maxx, maxy]])


class ParaxialModel():
    def __init__(self, opt_model):
        self.opt_model = opt_model
        self.seq_model = opt_model.seq_model
        self.sys = None
        self.ax = []
        self.pr = []
        self.opt_inv = 1.0

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']
        del attrs['seq_model']
        del attrs['parax_data']
        return attrs

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        self.seq_model = opt_model.seq_model

    def update_model(self):
        self.parax_data = self.opt_model.optical_spec.parax_data
        self.build_lens()

    def build_lens(self):
        # rebuild the `sys` description from the seq_model path
        self.sys = sys = self.seq_path_to_paraxial_lens(self.seq_model.path())

        # precalculate the reduced forms of the paraxial axial and chief rays
        ax_ray, pr_ray, fod = self.parax_data
        self.opt_inv = fod.opt_inv

        self.ax = []
        self.pr = []
        for i in range(0, len(sys)):
            n = sys[i][indx]
            self.ax.append([ax_ray[i][ht], n*ax_ray[i][slp], n*ax_ray[i][aoi]])
            self.pr.append([pr_ray[i][ht], n*pr_ray[i][slp], n*pr_ray[i][aoi]])

    def add_node(self, node, new_vertex, type_sel, interact_mode):
        """ Add a node in the paraxial data structures """
        ns = self.seq_model.get_num_surfaces()
        if node >= ns - 1:
            node = ns - 2
        n = self.sys[node][indx]
        new_node = node + 1
        self.sys.insert(new_node, [0.0, 0.0, n, interact_mode])

        ax_node = [0.0, 0.0, 0.0]
        ax_node[type_sel] = new_vertex[1]
        self.ax.insert(new_node, ax_node)

        pr_node = [0.0, 0.0, 0.0]
        pr_node[type_sel] = new_vertex[0]
        self.pr.insert(new_node, pr_node)

        if type_sel == ht:
            self.apply_ht_dgm_data(new_node, new_vertex=new_vertex)
        elif type_sel == slp:
            self.apply_slope_dgm_data(new_node, new_vertex=new_vertex)

        if interact_mode == 'reflect':
            for i in range(new_node, len(self.sys)):
                self.sys[i][indx] = -self.sys[i][indx]

        return new_node

    def assign_object_to_node(self, node, factory, **inputs):
        """ create a new element from `factory` and replace `node` with it """

        # extract optical properties of node
        n = self.sys[node][indx]
        power = self.sys[node][pwr]
        thi = n*self.sys[node][tau]
        sd = abs(self.ax[node][ht]) + abs(self.pr[node][ht])

        # create an element with the node's properties
        seq, ele = factory(power=power, sd=sd)

        n_before = self.sys[node-1][indx]
        thi_before = n_before*self.sys[node-1][tau]
        self.seq_model.gaps[node-1].thi = thi_before

        # insert the path sequence and elements into the
        #  sequential and element models
        args = seq, ele
        kwargs = {'idx': node-1, 't': thi, **inputs}
        self.opt_model.insert_ifc_gp_ele(*args, **kwargs)

        path_stop = node + len(seq)
        inserted_seq = list(self.seq_model.path(start=node-1, stop=path_stop))
        sys_seq = self.seq_path_to_paraxial_lens(inserted_seq[1:])
        pp_info = self.compute_principle_points(node, inserted_seq)
        self.replace_node_with_seq(node, sys_seq, pp_info)

        return args, kwargs

    def compute_principle_points(self, node, seq):
        n_0 = seq[0][mc.Indx]
        z_dir_before = seq[0][mc.Zdir]
        n_k = seq[-1][mc.Indx]
        z_dir_k = seq[-1][mc.Zdir]
        path = [[Surface(), Gap(), None, n_0, z_dir_before]]
        path.extend(seq[1:])
        pp_info = fo.compute_principle_points(iter(path),
                                              n_0=z_dir_before*n_0,
                                              n_k=z_dir_k*n_k)
        return pp_info

    def replace_node_with_seq(self, node, sys_seq, pp_info):
        """ replaces the data at node with sys_seq """
        efl, pp1, ppk, ffl, bfl = pp_info[2]
        self.sys[node-1][tau] -= pp1/self.sys[node-1][indx]
        sys_seq[-1][tau] = self.sys[node][tau] - ppk/sys_seq[-1][indx]

        self.delete_node(node)
        for ss in sys_seq:
            self.sys.insert(node, ss)
            self.ax.insert(node, [0.0, 0.0, 0.0])
            self.pr.insert(node, [0.0, 0.0, 0.0])
            node += 1
        self.paraxial_trace()

    def get_object_for_node(self, node):
        ''' basic 1:1 relationship between seq and parax model sequences '''
        ifc = self.seq_model.ifcs[node]
        ele = self.opt_model.ele_model.ifcs_dict[ifc]
        args = [[ifc, None, None, 1, 1]], [ele]
        kwargs = {'idx': node}
        return args, kwargs

    def add_object(self, surf, new_vertex, type_sel, factory, interact_mode):
        new_surf = self.add_node(surf, new_vertex, type_sel, interact_mode)
        self.assign_object_to_node(new_surf, factory)

    def delete_node(self, surf):
        """ delete the node at position surf """
        del self.sys[surf]
        del self.ax[surf]
        del self.pr[surf]

    def apply_ht_dgm_data(self, surf, new_vertex=None):
        """ This routine calculates all data dependent on the input
            height coordinates (y,ybar) at surface surf.
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr
        opt_inv = self.opt_inv

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

        if new_vertex is not None:
            pr_ray[surf][slp] = new_vertex[0]
            ax_ray[surf][slp] = new_vertex[1]

        nsm1 = len(sys) - 1
        if nsm1 == 0:
            p = 0
            c = 1
            ax_ray[c][ht] = ax_ray[p][slp]*sys[p][tau] + ax_ray[p][ht]
            pr_ray[c][ht] = pr_ray[p][slp]*sys[p][tau] + pr_ray[p][ht]

        else:
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

            if (surf < nsm1):
                s = surf + 1
                ax_ray[s][ht] = ax_ray[c][slp]*sys[c][tau] + ax_ray[c][ht]
                pr_ray[s][ht] = pr_ray[c][slp]*sys[c][tau] + pr_ray[c][ht]

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

    def list_lens(self):
        """ list the paraxial axial and chief rays, and power, reduced distance
        """
        sys = self.sys
        ax_ray = self.ax
        pr_ray = self.pr

        print("       ax_ray_ht    ax_ray_slp")
        for i in range(0, len(ax_ray)):
            print("{:2}: {:12.5g}  {:12.6g}".format(i, ax_ray[i][ht],
                                                    ax_ray[i][slp]))

        print("\n       pr_ray_ht    pr_ray_slp")
        for i in range(0, len(pr_ray)):
            print("{:2}: {:12.5g}  {:12.6g}".format(i, pr_ray[i][ht],
                                                    pr_ray[i][slp]))

        print("\n            power           tau        index    type")
        for i in range(0, len(sys)):
            print("{:2}: {:13.7g}  {:12.5g} {:12.5f}    {}".format(
                i, sys[i][pwr], sys[i][tau], sys[i][indx], sys[i][rmd]))

    def first_order_data(self):
        """List out the first order imaging properties of the model."""
        self.opt_model.optical_spec.parax_data.fod.list_first_order_data()

    def seq_path_to_paraxial_lens(self, path):
        """ returns lists of power, reduced thickness, signed index and refract
            mode
        """
        sys = []
        for i, sg in enumerate(path):
            n_after = sg[mc.Indx] if sg[mc.Zdir] > 0 else -sg[mc.Indx]
            imode = sg[mc.Intfc].interact_mode
            power = sg[mc.Intfc].optical_power
            if sg[mc.Gap]:
                tau = sg[mc.Gap].thi/n_after
                sys.append([power, tau, n_after, imode])
            else:
                sys.append([power, 0.0, n_after, imode])
        return sys

    def paraxial_lens_to_seq_model(self):
        """ Applies a paraxial lens spec (power, reduced distance) to the model
        """
        sys = self.sys
        ax_ray = self.ax
        n_before = sys[0][indx]
        slp_before = self.ax[0][slp]
        for i, sg in enumerate(self.seq_model.path()):
            if sg[mc.Gap]:
                n_after = sys[i][indx]
                slp_after = ax_ray[i][slp]
                sg[mc.Gap].thi = n_after*sys[i][tau]

                sg[mc.Intfc].set_optical_power(sys[i][pwr], n_before, n_after)
                sg[mc.Intfc].from_first_order(slp_before, slp_after,
                                              ax_ray[i][ht])

                n_before = n_after
                slp_before = slp_after

    def pwr_slope_solve(self, ray, surf, slp_new):
        p = ray[surf-1]
        c = ray[surf]
        pwr = (p[slp] - slp_new)/c[ht]
        return pwr

    def pwr_ht_solve(self, ray, surf, ht_new):
        sys = self.sys
        p = ray[surf-1]
        c = ray[surf]
        slp_new = (ht_new - c[ht])/sys[surf][tau]
        pwr = (p[slp] - slp_new)/ht_new
        return pwr

    def thi_ht_solve(self, ray, surf, ht_new):
        c = ray[surf]
        thi = (ht_new - c[ht])/c[slp]
        return thi

    def apply_data(self, vertex, lcl_pt):
        ray = self.ray
        p = ray[vertex-1]
        c = ray[vertex]
        c_slp_new = (lcl_pt[1] - c[ht])/lcl_pt[0]
        pwr = (p[slp] - c_slp_new)/c[ht]
        self.opt_model.seq_model.ifcs[vertex].optical_power = pwr
