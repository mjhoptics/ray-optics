#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright © 2021 Michael J. Hayford
"""Manage connectivity between sequence and element models using a tree.

.. Created on Sat Jan 23 20:15:48 2021

.. codeauthor: Michael J. Hayford
"""
from itertools import zip_longest

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

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        if hasattr(self, 'root_node'):
            root_node_compressed = self.root_node
            importer = DictImporter()
            self.root_node = importer.import_(root_node_compressed)
            sync_part_tree_on_restore(self.opt_model.ele_model,
                                      self.opt_model.seq_model,
                                      self.root_node)

    def update_model(self, **kwargs):
        sync_part_tree_on_update(self.opt_model.ele_model,
                                 self.opt_model.seq_model,
                                 self.root_node)
        self.sort_using_sequence(self.opt_model.seq_model)

    def init_from_sequence(self, seq_model):
        """Initialize part tree using a *seq_model*. """
        root_node = self.root_node
        for i, sgz in enumerate(zip_longest(seq_model.ifcs, seq_model.gaps,
                                            seq_model.z_dir)):
            s, gap, z_dir = sgz
            Node(f'i{i}', id=s, tag='#ifc', parent=root_node)
            if gap is not None:
                Node(f'g{i}', id=(gap, z_dir), tag='#gap', parent=root_node)

    def sort_using_sequence(self, seq_model):
        """Resequence part tree using a *seq_model*. """
        e_node_list = []
        for i, sgz in enumerate(zip_longest(seq_model.ifcs, seq_model.gaps,
                                            seq_model.z_dir)):
            ifc, gap, z_dir = sgz
            e_node = self.parent_node(ifc)
            if e_node is not None and e_node not in e_node_list:
                e_node_list.append(e_node)
            if gap is not None:
                g_node = self.parent_node((gap, z_dir))
                if g_node is not None and g_node not in e_node_list:
                    e_node_list.append(g_node)

        if len(e_node_list) > 0:
            self.root_node.children = e_node_list

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
        """ Return the node paired with `obj`. """
        return find_by_attr(self.root_node, name='id', value=obj)


    def trim_node(self, obj):
        """ Remove the branch where `obj` is the sole leaf. """
        leaf_node = self.node(obj)
        parent_node = None
        if leaf_node:
            parent_node = leaf_node.parent
        while parent_node is not None:
            if len(parent_node.children) > 1:
                # parent has more than one child, trim leaf_node
                leaf_node.parent = None
                break
            else:
                # trim leaf_node and continue up the branch
                leaf_node = parent_node
                parent_node = leaf_node.parent
                leaf_node.parent = None

    def parent_node(self, obj, tag='#element#airgap#dummyifc'):
        """ Return the parent node for `obj`, filtered by `tag`. """
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
        """ Return the parent object (and node) for `obj`, filtered by `tag`. """
        parent_node = self.parent_node(obj, tag)
        parent = parent_node.id if parent_node else None
        return parent, parent_node

    def list_tree(self, *args, **kwargs):
        """ Print a graphical console representation of the tree. 

        The optional arguments are passed through to the by_attr filter.
        Useful examples or arguments include:

            - pt.list_tree(lambda node: f"{node.name}: {node.tag}")
            - pt.list_tree(attrname='tag')

        """
        print(RenderTree(self.root_node).by_attr(*args, **kwargs))

    def list_tree_full(self):
        """ Print a graphical console representation of the tree with tags. """
        self.list_tree(lambda node: f"{node.name}: {node.tag}")

    def nodes_with_tag(self, tag='#element', not_tag='', root=None):
        """ Return a list of nodes that contain the requested `tag`. """
        def tag_filter(tags):
            def filter_tagged_node(node):
                for t in tags:
                    if t in node.tag:
                        for n_t in not_tags:
                            if n_t in node.tag:
                                return False
                        return True
                return False
            return filter_tagged_node

        tags = tag.split('#')[1:]
        not_tags = not_tag.split('#')[1:]
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
            try:
                idx = int(name[1:]) - 1
            except ValueError:
                idx = 0
            node.id = e.interface_list()[idx].profile
        elif name[:2] == 'tl':
            p_name = node.parent.name
            e = ele_dict[p_name]
            node.id = e.intrfc
        elif name[0] == 't':
            p_name = node.parent.name
            e = ele_dict[p_name]
            idx = int(name[1:])-1 if len(name) > 1 else 0
            node.id = e.gap_list()[idx]
        elif name[:1] == 'di':
            p_name = node.parent.name
            e = ele_dict[p_name]
            node.id = e.ref_ifc


def sync_part_tree_on_update(ele_model, seq_model, root_node):
    """Update node names to track element labels. """
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
        elif name == 'root':
            pass
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
        if part_tree.parent_node(ifc) is None:
            g_tfrm = g_tfrms[i]
    
            if g is not None:
                if g.medium.name().lower() == 'air':
                    num_eles = len(eles)
                    if num_eles == 0:
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
    
                            e_node = part_tree.add_element_to_tree(e, z_dir=z_dir1)
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
    
                            e_node = part_tree.add_element_to_tree(e, z_dir=z_dir)
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
                        ag_node = part_tree.add_element_to_tree(ag, z_dir=z_dir)
                        ag_node.leaves[0].id = (g, z_dir)
                        ele_model.add_element(ag)
    
                        eles = []
                        buried_reflector = False
    
                else:  # a non-air medium
                    # handle buried mirror, e.g. prism or Mangin mirror
                    if ifc.interact_mode == 'reflect':
                        buried_reflector = True
    
                    eles.append((i, ifc, g, z_dir, g_tfrm))
            else:
                process_airgap(ele_model, seq_model, part_tree,
                               i, g, z_dir, ifc, g_tfrm)

    # rename and tag the Image space airgap
    node = part_tree.parent_node((seq_model.gaps[-1], seq_model.z_dir[-1]))
    if node.name != 'Object space':
        node.name = node.id.label = 'Image space'
        node.tag += '#image'

    # sort the final tree by seq_model order
    part_tree.sort_using_sequence(seq_model)

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
    elif (s.interact_mode == 'dummy' or s.interact_mode == 'transmit'):
        add_dummy = False
        dummy_label = None
        dummy_tag = ''
        if i == 0:
            add_dummy = True  # add dummy for the object
            dummy_label = 'Object'
            dummy_tag = '#object'
        elif i == seq_model.get_num_surfaces() - 1:
            add_dummy = True  # add dummy for the object
            dummy_label = 'Image'
            dummy_tag = '#image'
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

    if g is not None:
        # add an AirGap
        if i == 0:
            gap_label = 'Object space'
            gap_tag = '#object'
        else:  # i > 0
            gap_label = None
            gap_tag = ''
        ag = elements.AirGap(g, idx=i, tfrm=g_tfrm, label=gap_label)
        ag_node = part_tree.add_element_to_tree(ag, z_dir=z_dir, tag=gap_tag)
        ag_node.leaves[0].id = (g, z_dir)
        ele_model.add_element(ag)
