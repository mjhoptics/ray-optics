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

import rayoptics.elem.elements as ele_module
import rayoptics.elem.sgz2ele as sgz2ele
import rayoptics.oprops.thinlens as thinlens
from rayoptics.util import str_to_class


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
        seq_model = self.opt_model.seq_model
        self.handle_object_image_tags(seq_model)
        self.sync_part_tree_on_update(self.opt_model.ele_model,
                                      seq_model,
                                      self.root_node)
        self.sort_tree_using_sequence(seq_model)

    def is_empty(self):
        if (isinstance(self.root_node, Node) and 
            len(self.root_node.children) == 0):
            return True
        else:
            return False

    def check_consistency(self, seq_model, ele_model):
        chg_list, sme, eme, ame = find_ele_changes(seq_model, ele_model, self)
        common_ele, added_ele, removed_ele, modified_ele = chg_list
        if len(added_ele)==0 and len(removed_ele)==0 and len(modified_ele)==0:
            return True
        else:
            return False

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
                    tag='#element#space#airgap#dummyifc#group')

                parse_path(path2root, groups)
        
            if gap is not None:
                gap_node = self.node((gap, z_dir))
                if gap_node is not None:
                    # path2root = gap_node.path[:-2]
                    path2root = self.nodes_with_tag(
                        node_list=gap_node.path, 
                        tag='#element#space#airgap#dummyifc#group')
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
    
    def remove_element_from_tree(self, e, **kwargs):
        e_node = self.node(e)
        e_node.children = []
        e_node.parent = None
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

    def parent_node(self, obj, tag='#element#space#dummyifc'):
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

    def parent_object(self, obj, tag='#element#space#dummyifc'):
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

    def build_pt_sg_lists(self):
        part_tag = '#assembly'
        nodes = self.nodes_with_tag(tag=part_tag)
        asms = [n.id for n in nodes]
        asm_list = []
        asm_dict = {}
        seq_model = self.opt_model['seq_model']
        for asm in asms:
            for p in asm.parts:
                ele_def = ele_module.build_ele_def(p, seq_model)
                asm_list.append(ele_def)
                asm_dict[ele_def] = asm
        return asm_list, asm_dict

    def list_pt_sg(self):
        ele_list, ele_dict = self.build_pt_sg_lists()
        for elem in ele_list:
            ele_type, idx_list, gap_list = elem
            e = ele_dict[elem]
            print(f"{e.label}: {ele_type[0]} {idx_list} {gap_list}")

    def sync_part_tree_on_update(self, ele_model, seq_model, root_node):
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

    def handle_object_image_tags(self, seq_model):
        """ Ensure nodes for object and image ifcs and gaps are tagged. """
        self._handle_object_image_tag(seq_model.ifcs[0], '#object', 
                                      not_tag='#space#airgap')
        self._handle_object_image_tag((seq_model.gaps[0], seq_model.z_dir[0]),
                                      '#object', not_tag='#dummyifc#surface')
        self._handle_object_image_tag(seq_model.ifcs[-1], '#image', 
                                      not_tag='#space#airgap')
        self._handle_object_image_tag((seq_model.gaps[-1],seq_model.z_dir[-1]),
                                      '#image', not_tag='#dummyifc#surface')

    def _handle_object_image_tag(self, sm_leaf_key, oi_tag, not_tag=''):
        nodes = self.nodes_with_tag(tag=oi_tag, not_tag=not_tag)
        oi_node = self.parent_node(sm_leaf_key)
        found_it = False
        for n in nodes:
            if n != oi_node:
                n.tag = n.tag.replace(oi_tag, '')
            else:
                found_it = True
        if not found_it and oi_node is not None:
            oi_node.tag += oi_tag


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


def sequence_to_elements(seq_model, ele_model, part_tree):
    """ Parse the seq_model into elements and update ele_model accordingly. """
    chg_list, sme, eme, ame = find_ele_changes(seq_model, ele_model, part_tree)
    common_ele, added_ele, removed_ele, modified_ele = chg_list
    sme_list, seq_str = sme
    eme_list, eme_dict = eme
    if len(added_ele)==0 and len(removed_ele)==0 and len(modified_ele)==0:
        # no additions, deletions or modifications to existing elements 
        # are necessary
        is_consistent = True

    else:
        # now the lists represent the 3 possible actions: adding, deleting,
        #  and modifying an existing element

        # elements in the added_ele list are created from the `ele_type`
        if len(added_ele) > 1:
            added_ele.sort(key=lambda ae: ae[1][0] 
                           if len(ae[1])>0 else ae[2][0])

        for ae in added_ele:
            (ele_token, ele_module, ele_class), *_ = ae
            e = str_to_class(ele_module, ele_class, 
                            ele_def_pkg=(seq_model, ae))
            ele_model.add_element(e)
            idx = e.reference_idx()
            z_dirs = seq_model.z_dir
            z_dir = z_dirs[idx] if idx < len(z_dirs) else z_dirs[idx-1]
            part_tree.add_element_to_tree(e, z_dir=z_dir)

        # items in the removed_ele list are removed from the ele_model 
        #  and part_tree
        for re in removed_ele:
            relem = eme_dict[re]
            ele_model.remove_element(relem)
            part_tree.remove_element_from_tree(relem)

        # modified elements use the `sync_to_ele_def` protocol to update an
        # existing element to a new ele_def
        for me in modified_ele:
            existing_ele, new_ele = me
            eme_dict[existing_ele].sync_to_ele_def(seq_model, new_ele)
        is_consistent = False

    return (is_consistent, chg_list, sme, eme)


def find_ele_changes(seq_model, ele_model, part_tree, print_visit=False):
    """ Parse the seq_model into elements and categorize the changes. 
    
    Returns:
        common_ele: list of ele_defs in common between sm and pt
        added_ele: list of ele_defs for new elements to be created 
        removed_ele: list of ele_defs to be removed 
        modified_ele: list of existing elements to be updated from new ele_defs  
        sme_list: ele_defs obtained by parsing the seq_model
        seq_str: character encoding of seq_model ifcs and gaps
        eme_list: ele_defs for current elements in the element model
        eme_dict: key: ele_def returns the value: element
        asm_list: ele_defs for current assemblies in the element/part model
        asm_dict: key: ele_def the value: assembly
    """
    # get sequential model "parse string"
    seq_str = seq_model.seq_str()
    if seq_str != '':
        sgz2ele_tree = sgz2ele.sgz2ele_grammar.parse(seq_str)
        
        sgz2ele_sm = sgz2ele.SMVisitor()
        sgz2ele_sm.do_print_visit = print_visit
        sgz2ele_visit = sgz2ele_sm.visit(sgz2ele_tree)
        
        sme_list = sgz2ele.flatten_visit(sgz2ele_visit)
        
        eme_list, eme_dict = ele_model.build_ele_sg_lists()
        
        eme_set = set(eme_list)
        sme_set = set(sme_list)
        
        common_ele = list(sme_set.intersection(eme_set))
        added_ele = list(sme_set.difference(eme_set))
        removed_ele = list(eme_set.difference(sme_set))

        # the modified element list is constructed from ele_defs with the same
        #  entity type and whose first gap indices match. 
        modified_ele = [(re, ae) for ae in added_ele for re in removed_ele 
                        if (ae[0][2] == re[0][2] and ae[2][0] == re[2][0])]
        # remove the modified_eles from their original lists 
        for me in modified_ele:
            re, ae = me
            added_ele.remove(ae)
            removed_ele.remove(re)    
        
        asm_list, asm_dict = part_tree.build_pt_sg_lists()

        return ((common_ele, added_ele, removed_ele, modified_ele), 
                (sme_list, seq_str),
                (eme_list, eme_dict),
                (asm_list, asm_dict))
    else:
        return (([], [], [], []), ([], ''), ([], {}), ([], {}))


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

