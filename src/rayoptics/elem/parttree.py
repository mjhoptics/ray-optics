#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2021 Michael J. Hayford
"""Manage connectivity between sequence and element models using a tree.

.. Created on Sat Jan 23 20:15:48 2021

.. codeauthor: Michael J. Hayford
"""

from anytree import Node, RenderTree, PreOrderIter
from anytree.exporter import DictExporter
from anytree.importer import DictImporter
from anytree.search import find_by_attr

from rayoptics.elem import elements
import rayoptics.oprops.thinlens as thinlens


class PartTree():
    def __init__(self, opt_model, **kwargs):
        self.opt_model = opt_model
        self.root_node = Node('root', id=self, tag='#group#root')

    def __json_encode__(self):
        attrs = dict(vars(self))
        del attrs['opt_model']

        exporter = DictExporter(attriter=lambda attrs_:
                                [(k, v) for k, v in attrs_
                                 if k == "name" or k == "tag"])
        attrs['root_node'] = exporter.export(self.root_node)
        return attrs

    def sync_to_restore(self):
        if hasattr(self, 'root_node'):
            root_node_compressed = self.root_node
            importer = DictImporter()
            self.root_node = importer.import_(root_node_compressed)
            sync_part_tree_on_restore(self.opt_model.ele_model,
                                      self.opt_model.seq_model,
                                      self.root_node)

    def update_model(self):
        sync_part_tree_on_update(self.opt_model.ele_model,
                                 self.opt_model.seq_model,
                                 self.root_node)

    def init_from_sequence(self, seq_model):
        """Initialize part tree using a *seq_model*. """
        for i, s in enumerate(seq_model.ifcs[1:-1], start=1):
            Node(f'i{i}', id=s, tag='#ifc', parent=self.root_node)
            gap = seq_model.gaps[i]
            z_dir = seq_model.z_dir[i]
            # if not isinstance(gap.medium, Air):
            Node(f'g{i}', id=(gap, z_dir), tag='#gap', parent=self.root_node)

    def add_element_model_to_tree(self, ele_model):
        for e in ele_model.elements:
            if hasattr(e, 'tree'):
                self.add_element_to_tree(e)
        return self

    def add_element_to_tree(self, e, **kwargs):
        e_node = e.tree(**kwargs)
        e_node.name = e.label
        leaves = e_node.leaves
        for leaf_node in leaves:
            dup_node = self.node(leaf_node.id)
            if dup_node is not None:
                dup_node.parent = None
        e_node.parent = self.root_node
        return e_node

    def node(self, obj):
        return find_by_attr(self.root_node, name='id', value=obj)

    def parent_node(self, obj, tag='#element#airgap#dummyifc'):
        tags = tag.split('#')[1:]
        leaf_node = self.node(obj)
        parent_node = leaf_node.parent if leaf_node else None
        while parent_node is not None:
            for t in tags:
                if t in parent_node.tag:
                    return parent_node
            parent_node = parent_node.parent
        return parent_node

    def parent_object(self, obj, tag='#element#airgap#dummyifc'):
        parent_node = self.parent_node(obj, tag)
        return parent_node.id if parent_node else None

    def list_tree(self):
        print(RenderTree(self.root_node).by_attr())

    def nodes_with_tag(self, tag='#element', root=None):
        def tag_filter(tags):
            def filter_tagged_node(node):
                for t in tags:
                    if t in node.tag:
                        return True
                return False
            return filter_tagged_node

        tags = tag.split('#')[1:]
        root_node = self.root_node if root is None else root
        nodes = [node for node in PreOrderIter(root_node,
                                               filter_=tag_filter(tags))]
        return nodes

    def list_model(self, tag='#element'):
        nodes = self.nodes_with_tag(tag)
        for i, node in enumerate(nodes):
            ele = node.id
            print("%d: %s (%s): %s" %
                  (i, ele.label, type(ele).__name__, ele))


def sync_part_tree_on_restore(ele_model, seq_model, root_node):
    ele_dict = {e.label: e for e in ele_model.elements}
    for node in PreOrderIter(root_node):
        name = node.name
        if name in ele_dict:
            node.id = ele_dict[name]
        elif name[0] == 'i':
            idx = int(name[1:])
            node.id = seq_model.ifcs[idx]
        elif name[0] == 'g':
            idx = int(name[1:])
            node.id = (seq_model.gaps[idx], seq_model.z_dir[idx])
        elif name[0] == 'p':
            p_name = node.parent.name
            e = ele_dict[p_name]
            idx = int(name[1:]) - 1
            node.id = e.interface_list()[idx].profile
        elif name[0] == 't':
            p_name = node.parent.name
            e = ele_dict[p_name]
            idx = int(name[1:]) - 1
            node.id = e.gap_list()[idx]
        elif name[:1] == 'di':
            p_name = node.parent.name
            e = ele_dict[p_name]
            node.id = e.ref_ifc
        elif name[:1] == 'tl':
            p_name = node.parent.name
            e = ele_dict[p_name]
            node.id = e.intrfc


def sync_part_tree_on_update(ele_model, seq_model, root_node):
    ele_dict = {e.label: e for e in ele_model.elements}
    for node in PreOrderIter(root_node):
        name = node.name
        if name[0] == 'i':
            idx = seq_model.ifcs.index(node.id)
            node.name = f'i{idx}'
        elif name[0] == 'g':
            gap, z_dir = node.id
            idx = seq_model.gaps.index(gap)
            z_dir = seq_model.z_dir[idx]
            node.id = (gap, z_dir)
            node.name = f'g{idx}'
        elif name[0] == 'p':
            p_name = node.parent.name
            e = ele_dict[p_name]
            idx = int(name[1:])-1 if len(name) > 1 else 0
            node.id = e.interface_list()[idx].profile
        elif name[:2] == 'di':
            p_name = node.parent.name
            e = ele_dict[p_name]
            node.id = e.ref_ifc
            idx = seq_model.ifcs.index(node.id)
            node.name = f'di{idx}'
        elif name[:2] == 'tl':
            p_name = node.parent.name
            e = ele_dict[p_name]
            node.id = e.intrfc
            idx = seq_model.ifcs.index(node.id)
            node.name = f'tl{idx}'
        elif name[0] == 't':
            p_name = node.parent.name
            e = ele_dict[p_name]
            idx = int(name[1:])-1 if len(name) > 1 else 0
            node.id = e.gap_list()[idx]
        else:
            if hasattr(node.id, 'label'):
                node.name = node.id.label


def elements_from_sequence(ele_model, seq_model, part_tree):
    """ generate an element list from a sequential model """

    if len(part_tree.root_node.children) == 0:
        # initialize part tree using the seq_model
        part_tree.init_from_sequence(seq_model)
    g_tfrms = seq_model.compute_global_coords(1)
    buried_reflector = False
    eles = []
    path = seq_model.path()
    for i, seg in enumerate(path):
        ifc, g, rindx, tfrm, z_dir = seg
        g_tfrm = g_tfrms[i]

        if g is not None:
            if g.medium.name().lower() == 'air':
                num_eles = len(eles)
                if num_eles == 0:
                    if i > 0:
                        process_airgap(
                            ele_model, seq_model, part_tree,
                            i, g, z_dir, ifc, g_tfrm, add_ele=True)
                else:
                    if buried_reflector is True:
                        num_eles = num_eles//2
                        eles.append((i, ifc, g, z_dir, g_tfrm))
                        i, ifc, g, z_dir, g_tfrm = eles[1]

                    if num_eles == 1:
                        i1, s1, g1, z_dir1, g_tfrm1 = eles[0]
                        sd = max(s1.surface_od(), ifc.surface_od())
                        e = elements.Element(s1, ifc, g1, sd=sd, tfrm=g_tfrm1,
                                             idx=i1, idx2=i)

                        e_node = part_tree.add_element_to_tree(e)
                        ele_model.add_element(e)
                        if buried_reflector is True:
                            ifc2 = eles[-1][1]
                            ifc_node = part_tree.node(ifc2)
                            p_node = find_by_attr(e_node, name='name',
                                                  value='p1')
                            ifc_node.parent = p_node

                            g1_node = part_tree.node((g1, z_dir1))
                            g_node = part_tree.node((g, z_dir))
                            g_node.parent = g1_node.parent

                            # set up for airgap
                            i, ifc, g, z_dir, g_tfrm = eles[-1]

                    elif num_eles > 1:
                        if not buried_reflector:
                            eles.append((i, ifc, g, z_dir, g_tfrm))
                        e = elements.CementedElement(eles[:num_eles+1])

                        e_node = part_tree.add_element_to_tree(e)
                        ele_model.add_element(e)
                        if buried_reflector is True:
                            for i, j in enumerate(range(-1, -num_eles-1, -1),
                                                  start=1):
                                ifc = eles[j][1]
                                ifc_node = part_tree.node(ifc)
                                pid = f'p{i}'
                                p_node = find_by_attr(e_node,
                                                      name='name',
                                                      value=pid)
                                ifc_node.parent = p_node
                                g = eles[j-1][2]
                                z_dir = eles[j-1][3]
                                g_node = part_tree.node((g, z_dir))
                                tid = f't{i}'
                                t_node = find_by_attr(e_node,
                                                      name='name',
                                                      value=tid)
                                if g_node:
                                    g_node.parent = t_node
                        # set up for airgap
                        i, ifc, g, z_dir, g_tfrm = eles[-1]

                    # add an AirGap
                    ag = elements.AirGap(g, idx=i, tfrm=g_tfrm)
                    ag_node = part_tree.add_element_to_tree(ag)
                    ag_node.leaves[0].id = (g, z_dir)
                    ele_model.add_element(ag)

                    eles = []
                    buried_reflector = False

            else:  # a non-air medium
                # handle buried mirror, e.g. prism or Mangin mirror
                if ifc.interact_mode == 'reflect':
                    buried_reflector = True

                eles.append((i, ifc, g, z_dir, g_tfrm))


def process_airgap(ele_model, seq_model, part_tree, i, g, z_dir, s, g_tfrm,
                   add_ele=True):
    if s.interact_mode == 'reflect' and add_ele:
        sd = s.surface_od()
        z_dir = seq_model.z_dir[i]
        m = elements.Mirror(s, sd=sd, tfrm=g_tfrm, idx=i, z_dir=z_dir)
        part_tree.add_element_to_tree(m)
        ele_model.add_element(m)
    elif isinstance(s, thinlens.ThinLens) and add_ele:
        te = elements.ThinElement(s, tfrm=g_tfrm, idx=i)
        part_tree.add_element_to_tree(te)
        ele_model.add_element(te)
    elif s.interact_mode == 'transmit':
        add_dummy = False
        dummy_label = None
        dummy_tag = ''
        if i == 0:
            add_dummy = True  # add dummy for the object
            dummy_label = 'Object'
            dummy_tag = '#object'
        else:  # i > 0
            gp = seq_model.gaps[i-1]
            if gp.medium.name().lower() == 'air':
                add_dummy = True
                if seq_model.stop_surface == i:
                    dummy_label = 'Stop'
                    dummy_tag = '#stop'
        if add_dummy:
            sd = s.surface_od()
            di = elements.DummyInterface(s, sd=sd, tfrm=g_tfrm, idx=i,
                                         label=dummy_label)
            part_tree.add_element_to_tree(di, tag=dummy_tag)
            ele_model.add_element(di)

    # add an AirGap
    ag = elements.AirGap(g, idx=i, tfrm=g_tfrm)
    ag_node = part_tree.add_element_to_tree(ag)
    ag_node.leaves[0].id = (g, z_dir)
    ele_model.add_element(ag)
