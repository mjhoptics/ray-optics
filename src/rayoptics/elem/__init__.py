""" Package providing support for element based optical modeling

    The :mod:`~.elem` subpackage provides classes and functions
    for element based models and rendering. These include:

        - Element based model, :mod:`~.elements`
        - Geometric constituents including :mod:`~.profiles` and
          :mod:`~.surface`
        - Lens layout and rendering support :mod:`~.layout`
        - Coordinate transformation support :mod:`~.transform`

    The element model is managed by the :class:`~.elements.ElementModel` class
"""

from anytree import Node

class RONode(Node):
    """ RayOptics subclass of Node 
    
    When deleting parts that are contained in Assembly objs, use the action of disconnecting the Part from the parttree to remove the Part from the assembly's list of parts.
    """
    def _pre_detach(self, parent):
        # print(f"{self.tag=},   {parent.tag=}")
        if "#assembly" in parent.tag:
            if self.id.parent is None:
                parent.id.parts.remove(self.id)
                # print(f"  Part {self.id.label} removed from assembly {parent.id.label}")
