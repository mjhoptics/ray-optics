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

    def trace(self, p0, d0, wl, eps=1.0e-12):
        path = nx.dfs_edges(self.sg)
        ray = []
        # trace object surface
        obj = next(path)
        srf = self.sg.node[obj[0]]['s']
        t, p = srf.profile.intersect(p0, d0)
        ray.append([p, d0])

        before = obj
        before_p = p
        before_d = d0

        # loop of remaining surfaces in path
        while True:
            try:
                after = next(path)
                before_gap = self.sg.get_edge_data(*before)['g']
                # t = before_gap.thi
                x = before_p - np.array([0, 0, before_gap.thi])
                # s = (before_p[2] - t)/before_d[2]
                eecd = -x.dot(before_d)
                # p_before = before_p + s*before_d
                before_eecp = x + eecd*before_d
                n_before = before_gap.medium.rindex(wl)
                srf = self.sg.node[before[1]]['s']
                after_gap = self.sg.get_edge_data(*after)['g']
                n_after = after_gap.medium.rindex(wl)
                # p, s = srf.profile.intersect(p_before, before_d)
                before_eecd, p = srf.profile.intersect(before_eecp, before_d)
                normal = srf.profile.normal(p)
                after_d = srf.profile.bend(before_d, normal, n_before, n_after)

                after_eecd = -p.dot(after_d)
                ray.append([p, after_d])

                before_p = p
                before_d = after_d
                before = after
                print(before_eecd, after_eecd)
            except StopIteration:
                break

        return ray
