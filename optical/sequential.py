import networkx as nx
import numpy as np
import globalspec
import surface as srf
import gap
import medium as m


class SequentialModel:
    def __init__(self):
        self.global_spec = globalspec.GlobalData()
        self.sg = nx.Graph()
        self.stop_surface = 1
        # self.cur_surface = 0
        self.cur_gap = (0, 1)
        self.next_surface = 2
        self.sg.add_node(self.cur_gap[0], s=srf.Surface('Obj'))
        self.sg.add_node(self.cur_gap[1], s=srf.Surface('Img'))
        self.sg.add_edge(*self.cur_gap, g=gap.Gap())

    def insert(self, surf, gap):
        """ insert surf and gap at the cur_gap edge of the sequential model
            graph """
        old_gap = self.sg.get_edge_data(*self.cur_gap)['g']
        self.sg.remove_edge(*self.cur_gap)
        self.sg.add_node(self.next_surface, s=surf)
        self.sg.add_edge(self.cur_gap[0], self.next_surface, g=old_gap)
        self.sg.add_edge(self.next_surface, self.cur_gap[1], g=gap)
        self.cur_gap = (self.next_surface, self.cur_gap[1])
        self.next_surface += 1

    def add_surface(self, surf):
        """ add a surface where surf is a list that contains:
            [curvature, thickness, refractive_index, v-number] """
        s = srf.Surface()
        s.profile.cv = surf[0]

        if surf[3] == 0.0:
            if surf[2] == 1.0:
                mat = m.Air()
            else:
                mat = m.Medium(surf[2])
        else:
            mat = m.Glass(surf[2], surf[3], '')

        g = gap.Gap(surf[1], mat)
        self.insert(s, g)

    def list_model(self):
        for gp in list(nx.dfs_edges(self.sg)):
            print(self.sg.node[gp[0]]['s'])
            print('   ', self.sg.get_edge_data(*gp)['g'])
        print(self.sg.node[gp[1]]['s'])

    def list_gaps(self):
        for gp in list(nx.dfs_edges(self.sg)):
            print(self.sg.get_edge_data(*gp)['g'])

    def list_surfaces(self):
        for gp in list(nx.dfs_edges(self.sg)):
            print(self.sg.node[gp[0]]['s'])
        print(self.sg.node[gp[1]]['s'])

    def trace(self, pt0, dir0, wl, eps=1.0e-12):
        path = nx.dfs_edges(self.sg)
        ray = []
        # trace object surface
        obj = next(path)
        srf_obj = self.sg.node[obj[0]]['s']
        gap_obj = self.sg.get_edge_data(*obj)['g']
        op_delta, pt = srf_obj.profile.intersect(pt0, dir0)
        ray.append([pt, dir0, op_delta])
        op_delta *= gap_obj.medium.rindex(wl)

        before = obj
        before_pt = pt
        before_dir = dir0

        # loop of remaining surfaces in path
        while True:
            try:
                after = next(path)

                gap_before = self.sg.get_edge_data(*before)['g']
                b4_pt_cur_coord = before_pt - np.array([0, 0, gap_before.thi])

                eic_dst = -b4_pt_cur_coord.dot(before_dir)

                eic_pt_before = b4_pt_cur_coord + eic_dst*before_dir

                n_before = gap_before.medium.rindex(wl)
                srf = self.sg.node[before[1]]['s']
                gap_after = self.sg.get_edge_data(*after)['g']
                n_after = gap_after.medium.rindex(wl)

                eic_dst_before, pt = srf.profile.intersect(eic_pt_before,
                                                           before_dir)
                normal = srf.profile.normal(pt)
                after_dir = srf.profile.bend(before_dir, normal,
                                             n_before, n_after)

                eic_dst_after = -pt.dot(after_dir)
                op_delta += (n_after*eic_dst_after - n_before*eic_dst_before)
                dst_before = eic_dst + eic_dst_before
                ray.append([pt, after_dir, dst_before])

                before_pt = pt
                before_dir = after_dir
                before = after
                print(eic_dst_before, eic_dst_after)
            except StopIteration:
                break

        return ray, op_delta
