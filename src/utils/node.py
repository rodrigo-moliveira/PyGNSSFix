""" Node Module
This module helps to find the shortest path between two elements in a hierarchy or in a graph.

This module is imported and adapted from the `Beyond Python Package <https://pypi.org/project/beyond/>`__.

See the original implementation `here <https://github.com/galactics/beyond/blob/master/beyond/utils/node.py>`__.

All credits to Jules David, the owner of Beyond library.
"""

from collections import OrderedDict


class Route:
    """
    Class used by :py:class:`Node` to describe where to find
    another node.
    """

    def __init__(self, direction, steps):
        self.direction = direction
        self.steps = steps

    def __repr__(self):  # pragma: no cover
        return f"<d={self.direction}, s={self.steps}>"


class Node:
    """
    Class representing a node in a graph, relations may be circular.

    Attributes:
        name(str): Name of the node
        neighbors(OrderedDict): List of all direct neighbors in the graph.
            OrderedDict is only used as OrderedSet, so only the keys of the dict matter
        routes(dict): Route mapping. What direction to follow in order to reach a particular target
    """

    def __init__(self, name):
        """
        Args:
            name (str): Name of the node. Will be used for graph searching
        """
        self.name = name

        self.neighbors = OrderedDict()

        self.routes = {}

    def __add__(self, other):
        self.neighbors[other] = None
        other.neighbors[self] = None
        self._update()
        return other

    @property
    def list(self):
        return [self.path(node_name)[-1] for node_name in self.routes.keys()] + [self]

    def _update(self, already_updated=None):

        self.routes = {}
        for node in self.neighbors:
            self.routes[node.name] = Route(node, 1)

            # Retrieve route from neighbors
            for name, route in node.routes.items():

                # check if the node actually at hand (name) is not already
                # a direct neighbor of self
                if name in [self.name] + [x.name for x in self.neighbors]:
                    continue

                # check if the node actually at hand (name) is not already
                # integrated in the self.routes or if it already is, if the
                # path is shorter
                if (
                    name in self.routes.keys()
                    and self.routes[name].steps <= route.steps
                ):
                    continue

                self.routes[name] = Route(node, route.steps + 1)

        # This set serves as a shared lock, every object that is in this set
        # won't be updated by others. This is to avoid infinite recursions
        if already_updated is None:
            already_updated = set()

        already_updated.add(self)

        # Recursive update (with lock)
        for node in self.neighbors:
            if node not in already_updated:
                node._update(already_updated)

    def path(self, goal):
        """
        Get the shortest way between two nodes of the graph.

        Args:
            goal (str): Name of the targeted node
        Returns:
            list[Node]: a list of the Nodes of the path
        Raises:
            ValueError: an exception is thrown if the goal node is unknown
        """

        if isinstance(goal, Node):
            goal = goal.name

        if goal == self.name:
            return [self]

        if goal not in self.routes:
            raise ValueError(f"Unknown '{goal}'")

        obj = self
        path = [obj]
        while True:
            obj = obj.routes[goal].direction
            path.append(obj)
            if obj.name == goal:
                break
        return path

    def steps(self, goal):
        """
        Get the list of individual relations leading to the targeted node.

        Args:
            goal (str): Name of the targeted node
        Return:
            list[Node]: a list of the Nodes of the path
        """

        path = self.path(goal)
        for i in range(len(path) - 1):
            yield path[i], path[i + 1]

    def __str__(self):  # pragma: no cover
        return self.name

    def __repr__(self):  # pragma: no cover
        return f"<{self.__class__.__name__} '{self.name}' at '{hex(id(self))}'>"
