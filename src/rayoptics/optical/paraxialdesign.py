#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2018 Michael J. Hayford
""" First order paraxial design space

.. Created on Sat Mar 31 21:14:42 2018

.. codeauthor: Michael J. Hayford
"""

from collections import namedtuple
import numpy as np

from rayoptics.optical.model_constants import ht, slp, aoi
from rayoptics.optical.model_constants import pwr, tau, indx, rmd
import rayoptics.optical.model_constants as mc
from rayoptics.optical.elements import insert_ifc_gp_ele
import rayoptics.optical.firstorder as fo
from rayoptics.optical.gap import Gap
from rayoptics.optical.surface import Surface
from rayoptics.optical.surface import InteractionMode as imode
from rayoptics.util.rgb2mpl import rgb2mpl, backgrnd_color


def bbox_from_poly(poly):
    minx, miny = np.min(poly, axis=0)
    maxx, maxy = np.max(poly, axis=0)
    return np.array([[minx, miny], [maxx, maxy]])


ParaxData = namedtuple('ParaxData', ['ht', 'slp', 'aoi'])
""" paraxial ray data at an interface

    Attributes:
        ht: height at interface
        slp: n*slope
        aoi: n*angle of incidence
"""

ParaxSys = namedtuple('ParaxSys', ['pwr', 'tau', 'indx', 'rmd'])
""" paraxial lens data at an interface

    Attributes:
        pwr: power
        tau: reduced distance, t/n
        indx: refractive index (n)
        rmd: refract mode
"""


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

    def update_model(self):
        self.parax_data = self.opt_model.optical_spec.parax_data
        self.build_lens()

    def build_lens(self):
        self.sys = self.seq_path_to_paraxial_lens(self.seq_model.path())
        sys = self.sys
        ax_ray, pr_ray, fod = self.parax_data
        self.opt_inv = fod.opt_inv

        self.ax = []
        self.pr = []
        for i in range(0, len(sys)):
            n = sys[i][indx]
            self.ax.append([ax_ray[i][ht], n*ax_ray[i][slp], n*ax_ray[i][aoi]])
            self.pr.append([pr_ray[i][ht], n*pr_ray[i][slp], n*pr_ray[i][aoi]])

    def add_node(self, surf, new_vertex, type_sel):
        """ Add a node in the paraxial data structures """
        ns = self.seq_model.get_num_surfaces()
        if surf >= ns - 1:
            surf = ns - 2
        n = self.sys[surf][indx]
        new_surf = surf + 1
        self.sys.insert(new_surf, [0.0, 0.0, n, imode.Transmit])

        ax_node = [0.0, 0.0, 0.0]
        ax_node[type_sel] = new_vertex[1]
        self.ax.insert(new_surf, ax_node)

        pr_node = [0.0, 0.0, 0.0]
        pr_node[type_sel] = new_vertex[0]
        self.pr.insert(new_surf, pr_node)

        if type_sel == ht:
            self.apply_ht_dgm_data(new_surf, new_vertex=new_vertex)
        elif type_sel == slp:
            self.apply_slope_dgm_data(new_surf, new_vertex=new_vertex)
        return new_surf

    def assign_object_to_node(self, node, factory):
        """ create a new element from `factory` and replace `node` with it """

        # extract optical properties of node
        n = self.sys[node][indx]
        power = self.sys[node][pwr]
        thi = n*self.sys[node][tau]
        sd = abs(self.ax[node][ht]) + abs(self.pr[node][ht])

        # create an element with the node's properties
        seq, ele = factory(power=power, sd=sd)
        # insert the path sequence and elements into the
        #  sequential and element models
        insert_ifc_gp_ele(self.opt_model, seq, ele, idx=node-1, t=thi)

        self.sys[node][rmd] = seq[0][0].interact_mode
        if seq[0][0].interact_mode == imode.Reflect:
            self.sys[node][indx] = -self.sys[node][indx]

        self.replace_node_with_seq(node, seq)

    def replace_node_with_seq(self, node, seq):
        """ replaces the data at node with seq """
        n_0 = self.sys[node-1][indx]
        z_dir_before = 1 if n_0 > 0 else -1
        n_k = seq[-1][mc.Indx]
        path = [[Surface(), Gap(), None, n_0, z_dir_before]]
        path.extend(seq)
        pp_info = fo.compute_principle_points(iter(path), n_k=n_k)
        efl, pp1, ppk, ffl, bfl = pp_info[2]
        seq_sys = self.seq_path_to_paraxial_lens(iter(seq))
        self.sys[node-1][tau] -= pp1/n_0
        seq_sys[-1][tau] = self.sys[node][tau] - ppk/seq_sys[-1][indx]
        self.delete_node(node)
        for ss in seq_sys:
            self.sys.insert(node, ss)
            self.ax.insert(node, [0.0, 0.0, 0.0])
            self.pr.insert(node, [0.0, 0.0, 0.0])
            node += 1
        self.paraxial_trace()

    def add_object(self, surf, new_vertex, type_sel, factory):
        new_surf = self.add_node(surf, new_vertex, type_sel)
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
            print("{:2}: {:13.7g}  {:12.5g} {:12.5f}    {}".format(i,
                  sys[i][pwr], sys[i][tau], sys[i][indx], sys[i][rmd].name))

    def seq_path_to_paraxial_lens(self, path):
        """ returns lists of power, reduced thickness, signed index and refract
            mode
        """
        sys = []
        for i, sg in enumerate(path):
            n_after = sg[mc.Indx] if sg[mc.Zdir] > 0 else -sg[mc.Indx]
            rmode = sg[mc.Intfc].interact_mode
            power = sg[mc.Intfc].optical_power
            if sg[mc.Gap]:
                tau = sg[mc.Gap].thi/n_after
                sys.append([power, tau, n_after, rmode])
            else:
                sys.append([power, 0.0, n_after, rmode])
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


class Diagram():
    """ class for paraxial ray rendering/editing """
    def __init__(self, opt_model, dgm_type, seq_start=1,
                 label='paraxial'):
        self.label = label
        self.opt_model = opt_model

        self.setup_dgm_type(dgm_type)

        if opt_model.specsheet.conjugate_type == 'finite':
            self.seq_start = 0
        elif opt_model.specsheet.conjugate_type == 'infinite':
            self.seq_start = 1
        else:
            self.seq_start = seq_start

        self.handles = {}
        self.actions = self.edit_diagram_actions()
        self.cur_node = None

    def setup_dgm_type(self, dgm_type):
        parax_model = self.opt_model.parax_model
        if dgm_type == 'ht':
            self.type_sel = ht
            self.apply_data = parax_model.apply_ht_dgm_data
        elif dgm_type == 'slp':
            self.type_sel = slp
            self.apply_data = parax_model.apply_slope_dgm_data

    def get_label(self):
        return self.label

    def register_commands(self, *args, **kwargs):
        fig = kwargs.pop('figure')
        cmd_actions = kwargs.pop('cmd_actions')
        actions = cmd_actions(**kwargs)

        def do_command_action(event, target, event_key):
            nonlocal fig, actions
            if target is not None:
                shape, handle = target.artist.shape
                try:
                    action = actions[event_key]
                    action(fig, handle, event, target.info)
                except KeyError:
                    pass
        fig.do_action = do_command_action

    def render_handles(self):
        parax_model = self.opt_model.parax_model
        self.handles = {}

        shape = []
        for x, y in zip(parax_model.pr, parax_model.ax):
            shape.append([x[self.type_sel], y[self.type_sel]])

        n_color = rgb2mpl([138, 43, 226])  # blueviolet
        self.handles['nodes'] = shape, 'polyline', {'linestyle': '',
                                                    'marker': 's',
                                                    'picker': 6,
                                                    'color': n_color,
                                                    'hilite': n_color}
        e_color = rgb2mpl([198, 113, 113])  # sgi salmon
        self.handles['edges'] = shape, 'polyline', {'color': e_color,
                                                    'hilite': e_color}

        return self.handles

    def handle_actions(self):
        self.actions = {}

        node_actions = {}
#        node_actions['pt'] = BendAction(self)
        self.actions['nodes'] = node_actions

        edge_actions = {}
#        edge_actions['pt'] = SagAction(self.s1)
        self.actions['edges'] = edge_actions

        return self.actions

    def update_shape(self, view):
        def gen_poly_artist(poly, key, kwargs):
            rp = view.create_polyline(poly, linewidth=1,
                                      **kwargs)
            return rp, bbox_from_poly(poly)

        handles = self.render_handles()

        gui_handles = {}
        for key, graphics_handle in handles.items():
            poly = np.array(graphics_handle[0])
            if graphics_handle[1] == 'polyline':
                gui_handles[key] = gen_poly_artist(poly, key,
                                                   graphics_handle[2])

        return gui_handles

    def edit_diagram_actions(self):
        parax_model = self.opt_model.parax_model
        actions = {}

        def on_select_point(fig, handle, event, info):
            self.cur_node = self.seq_start
            if 'ind' in info:
                self.cur_node += info['ind'][0]
        actions['press'] = on_select_point

        def on_edit_point(fig, handle, event, info):
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        actions['drag'] = on_edit_point

        def on_release_point(fig, handle, event, info):
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
            self.cur_node = None
        actions['release'] = on_release_point

        return actions

    def add_element_actions(self, factory=None, node_init=None, **kwargs):
        parax_model = self.opt_model.parax_model
        actions = {}

        def on_press_add_point(fig, handle, event, info):
            self.cur_node = self.seq_start
            if 'ind' in info:
                self.cur_node += info['ind'][0]
            print("add_point: press", self.cur_node)
            event_data = np.array([event.xdata, event.ydata])
            parax_model.add_node(self.cur_node, event_data, self.type_sel)
            self.cur_node += 1
            parax_model.assign_object_to_node(self.cur_node, node_init)
            fig.parax_model.paraxial_lens_to_seq_model()
            fig.skip_build = True
            fig.refresh_gui()
        actions['press'] = on_press_add_point

        def on_drag_add_point(fig, handle, event, info):
            print("add_point: drag", self.cur_node)
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        actions['drag'] = on_drag_add_point

        def on_release_add_point(fig, handle, event, info):
            print("add_point: release", self.cur_node)
            parax_model.assign_object_to_node(self.cur_node, factory)
            parax_model.paraxial_lens_to_seq_model()
            fig.skip_build = True
            fig.refresh_gui()
            self.cur_node = None
        actions['release'] = on_release_add_point

        return actions



class EditNodeAction():
    """ Action to move a diagram node, using an input pt """
    def __init__(self, dgm, parax_model, dgm_type):
        self.parax_model = parax_model
        self.select_pt_lcl = None
        self.cur_value = None
        self.new_value = None
        self.actions = {}

        def on_select(fig, event):
            self.cur_node = self.seq_start
            if 'ind' in info:
                self.cur_node += info['ind'][0]
        self.actions['press'] = on_select

        def on_edit(fig, event, pt):
            event_data = np.array([event.xdata, event.ydata])
            self.apply_data(self.cur_node, event_data)
            self.parax_model.paraxial_lens_to_seq_model()
            fig.refresh_gui()
        self.actions['drag'] = on_edit

        def on_release(fig, event):
            self.cv_new = self.ele.reference_interface().profile_cv
            xsag = sag(self.cv_new, event.lcl_pt[1])
#            print('on_release: {:.3f} {:.3f} {:.5f}'
#                  .format(event.lcl_pt[0], xsag, self.cv_new))
            fig.refresh_gui()
        self.actions['release'] = on_release
