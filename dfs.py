"""Adapted from NetworkX DFS on Edges (NOT written by me)"""

FORWARD = 'forward'
REVERSE = 'reverse'


def helper_funcs(G, orientation):
    """
    These are various G-specific functions that help us implement the algorithm
    for all graph types: graph, multigraph, directed or not.

    """
    ignore_orientation = G.is_directed() and orientation == 'ignore'
    reverse_orientation = G.is_directed() and orientation == 'reverse'

    if ignore_orientation:
        # When we ignore the orientation, we still need to know how the edge
        # was traversed, so we add an object representing the direction.
        def out_edges(u, **kwds):
            for edge in G.out_edges(u, **kwds):
                yield edge + (FORWARD,)
            for edge in G.in_edges(u, **kwds):
                yield edge + (REVERSE,)
    elif reverse_orientation:
        def out_edges(u, **kwds):
            for edge in G.in_edges(u, **kwds):
                yield edge + (REVERSE,)
    else:
        # If "yield from" were an option, we could pass kwds automatically.
        out_edges = G.edges

    # If every edge had a unique key, then it would be easier to track which
    # edges had been visited. Since that is not available, we will form a
    # unique identifier from the edge and key (if present). If the graph
    # is undirected, then the head and tail need to be stored as a frozenset.
    if ignore_orientation or reverse_orientation:
        # edge is a 4-tuple: (u, v, key, direction)
        # u and v always represent the true tail and head of the edge.
        def key(edge):
            # We want everything but the direction.
            return edge[:-1]
    else:
        if G.is_directed():
            def key(edge):
                return edge
        else:
            # edge is a 3-tuple:  (u, v, key)
            def key(edge):
                new_edge = (frozenset(edge[:2]),) + edge[2:]
                return new_edge

    def traversed_tailhead(edge):
        """
        Returns the tail and head of an edge, as it was traversed.

        So in general, this is different from the true tail and head.
        (Also, undirected edges have no true tail or head.)

        """
        if (ignore_orientation or reverse_orientation) and edge[-1] == REVERSE:
            tail, head = edge[1], edge[0]
        else:
            tail, head = edge[0], edge[1]
        return tail, head

    return out_edges, key, traversed_tailhead

def custom_edge_dfs(G, start_node, orientation='original'):
    """
    Depth first search
    """
    kwds = {'data': False}
    if G.is_multigraph():
        kwds['keys'] = True

    out_edges, key, tailhead = helper_funcs(G, orientation)
    
    goal_path_length = len(G.edges)
    stack = [start_node]
    node_edges = {}
    tried_paths = set()
    path = []
    paths = []
    
    while stack:

        # Get a complete list of all edges from the current node
        current_node = stack[-1]
        if current_node not in node_edges:
            node_edges[current_node] = list(out_edges(current_node, **kwds))
        
        # Enumerate the edges that are not currently used and have not been
        # tried
        check_edges = []
        for edge in node_edges[current_node]:
            if edge in path:
                continue
            chk = tuple(path + [edge])
            if chk not in tried_paths:
                check_edges.append(edge)
        
        # Try next edge, if available
        try:
            edge = check_edges.pop()
            stack.append(edge[1])
            path.append(edge)
            tried_paths.add(tuple(path))
        except IndexError:
            stack.pop()
            if path:
                rmedge = path.pop()
        
        # If we have used all edges then we're done
        if len(path) == goal_path_length:
            paths.append(list(path))
    
    return paths
