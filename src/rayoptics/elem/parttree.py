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

        def node_attrs(attrs_):
            attrs = []
            for k, v in attrs_:
                if k == "name" or k == "tag":
                    attrs.append((k, v))
                elif k == "id":
                    attrs.append(("id_key", str(id(v))))
            return attrs

        exporter = DictExporter(attriter=node_attrs)
        attrs['root_node'] = exporter.export(self.root_node)
        return attrs

    def sync_to_restore(self, opt_model):
        self.opt_model = opt_model
        if hasattr(self, 'root_node'):
            root_node_compressed = self.root_node
            importer = DictImporter()
            self.root_node = importer.import_(root_node_compressed)
            self.root_node.id = self
            if hasattr(self.root_node, 'id_key'):
                sync_part_tree_on_restore_idkey(self.opt_model, 
                                                self.opt_model.ele_model,
                                                self.opt_model.seq_model,
                                                self.root_node)
            else:
                sync_part_tree_on_restore(self.opt_model, 
                                          self.opt_model.ele_model,
                                          self.opt_model.seq_model,
                                          self.root_node)

    def update_model(self, **kwargs):
        sync_part_tree_on_update(self.opt_model.ele_model,
                                 self.opt_model.seq_model,
                                 self.root_node)
        self.sort_tree_using_sequence(self.opt_model.seq_model)

    def init_from_sequence(self, seq_model):
        """Initialize part tree using a *seq_model*. """
        root_node = self.root_node
        for i, sgz in enumerate(zip_longest(seq_model.ifcs, seq_model.gaps,
                                            seq_model.z_dir)):
            s, gap, z_dir = sgz
            Node(f'i{i}', id=s, tag='#ifc', parent=root_node)
            if gap is not None:
                Node(f'g{i}', id=(gap, z_dir), tag='#gap', parent=root_node)

    def sort_tree_using_sequence(self, seq_model):
        """Resequence part tree using a *seq_model*. """
        def parse_path(path2root, groups):
            for i, node in enumerate(path2root[1:]):
                parent_node = path2root[i]
                if parent_node in groups:
                    if node not in groups[parent_node]:
                        groups[parent_node].append(node)
                else:
                    groups[parent_node] = [node]

        groups = {}
        for i, sgz in enumerate(zip_longest(seq_model.ifcs, 
                                            seq_model.gaps, 
                                            seq_model.z_dir)):        
            ifc, gap, z_dir = sgz
        
            ifc_node = self.node(ifc)
            if ifc_node is not None:
                # path2root = ifc_node.path[:-2]
                path2root = self.nodes_with_tag(
                    node_list=ifc_node.path, 
                    tag='#element#airgap#dummyifc#group')

                parse_path(path2root, groups)
        
            if gap is not None:
                gap_node = self.node((gap, z_dir))
                if gap_node is not None:
                    # path2root = gap_node.path[:-2]
                    path2root = self.nodes_with_tag(
                        node_list=gap_node.path, 
                        tag='#element#airgap#dummyifc#group')
                    parse_path(path2root, groups)

        for group_node, child_list in groups.items():
            group_node.children = child_list

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

    def obj_by_name(self, name):
        """ Return the node paired with `obj`. """
        return find_by_attr(self.root_node, name='name', value=name).id

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

    def get_child_filter(self, tag='#element#assembly', not_tag=''):
        """ Returns a fct that filters a list of nodes to satisfy the tags"""
        def children_with_tag(children):
            return self.nodes_with_tag(node_list=children, 
                                       tag=tag, not_tag=not_tag)
        return children_with_tag

    def list_tree(self, *args, **kwargs):
        """ Print a graphical console representation of the tree. 

        The optional arguments are passed through to the by_attr filter.
        Useful examples or arguments include:

            - pt.list_tree(lambda node: f"{node.name}: {node.tag}")
            - pt.list_tree(attrname='tag')

        """
        list_tree_from_node(self.root_node, *args, **kwargs)

    def list_tree_full(self):
        """ Print a graphical console representation of the tree with tags. """
        self.list_tree(lambda node: f"{node.name}: {node.tag}")

    def nodes_with_tag(self, tag='#element', not_tag='',
                       root=None, node_list=None):
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
        
        if root is not None:
            nodes = [node for node in PreOrderIter(root,
                                                   filter_=tag_filter(tags))]
        elif node_list is not None:
            filter_node = tag_filter(tags)
            nodes = [node for node in node_list if filter_node(node)]
        else:
            root_node = self.root_node
            nodes = [node for node in PreOrderIter(root_node,
                                                   filter_=tag_filter(tags))]
        return nodes

    def list_model(self, tag='#element#assembly#dummyifc'):
        self.list_tree(childiter=self.get_child_filter(tag=tag))


def sync_part_tree_on_restore(opt_model, ele_model, seq_model, root_node):
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
            node.id = e.profile_list()[idx]
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


def sync_part_tree_on_restore_idkey(opt_model, ele_model, seq_model, root_node):
    for node in PreOrderIter(root_node):
        name = node.name
        if node.id_key in opt_model.parts_dict:
            node.id = opt_model.parts_dict[node.id_key]
        elif name[0] == 'i':
            idx = int(name[1:])
            node.id = seq_model.ifcs[idx]
        elif name[0] == 'g':
            idx = int(name[1:])
            node.id = (seq_model.gaps[idx], seq_model.z_dir[idx])
        elif name[0] == 'p':
            node.id = opt_model.profile_dict[node.id_key]
        elif name[:2] == 'tl':
            e = opt_model.parts_dict[node.parent.id_key]
            node.id = e.intrfc
        elif name[0] == 't':
            e = opt_model.parts_dict[node.parent.id_key]
            idx = int(name[1:])-1 if len(name) > 1 else 0
            node.id = e.gap_list()[idx]
        elif name[:1] == 'di':
            e = opt_model.parts_dict[node.parent.id_key]
            node.id = e.ref_ifc

    for node in PreOrderIter(root_node):
        delattr(node, 'id_key')


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
            num_idxs = len(e.idx_list())
            idx = (num_idxs-idx-1) if e.is_flipped else idx
            node.id = e.profile_list()[idx]
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
            if hasattr(node, 'id'):
                if hasattr(node.id, 'label'):
                    node.name = node.id.label
            else:
                print(f"sync_part_tree_on_update: No id attribute: {node.name}, {node.tag}")


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
        ifc, g, tfrm, rindx, z_dir = seg
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
    part_tree.sort_tree_using_sequence(seq_model)


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


def part_list_from_seq(opt_model, idx1, idx2):
    """Using the part_tree, return the parts for the input sequence range. """
    sm = opt_model['seq_model']
    seq = zip_longest(sm.ifcs[idx1:idx2+1], 
                      sm.gaps[idx1:idx2], 
                      sm.z_dir[idx1:idx2])
    part_tree = opt_model['part_tree']
    node_set = set()
    for ifc, gap, z_dir in seq:
        ifc_node = part_tree.node(ifc)
        parent_asm = ifc_node.ancestors[1]
        node_set.add(parent_asm)
        if gap is not None:
            gap_node = part_tree.node((gap, z_dir))
            parent_asm = gap_node.ancestors[1]
            node_set.add(parent_asm)

    node_list = sorted(node_set, key=lambda node: node.id.reference_idx())
    part_list = [node.id for node in node_list]
    return part_list, node_list


def list_tree_all_from_node(node, **kwargs):
    """ List the tree from `node` with full node output. """
    tag_filter = kwargs.pop('childiter', list)
    print(RenderTree(node, childiter=tag_filter))


def list_tree_from_node(node, *args, **kwargs):
    """ List the tree from `node` with attribute filtering. 

    The optional arguments are passed through to the by_attr filter.
    Useful examples or arguments include:

        - pt.list_tree(lambda node: f"{node.name}: {node.tag}")
        - pt.list_tree(attrname='tag')

    """
    tag_filter = kwargs.pop('childiter', list)
    print(RenderTree(node, childiter=tag_filter).by_attr(*args, **kwargs))
