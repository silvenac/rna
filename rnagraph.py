import re
import networkx as nx
from digest import (to_digest, to_graph, interior_bases)
from dfs import custom_edge_dfs
from networkx.drawing.nx_agraph import write_dot


class RNAGraph(object):
    """ A class for RNA Graph: instantiates with g and uc fragmentations of an
    RNA sequence.

    self.graph: NetworkX MultiDigraph used to calculate Eulerian Trails
    self.eulergraph: NetworkX MultiDigraph Eulerian representation of RNA Seq
    self.start: Start point of RNA sequence
    self.end: End point of RNA sequence

    In PyGraphViz graph drawing, start node is indicated by + and end node by *
    """

    def __init__(self, g_frag, uc_frag):
        self.start = None
        self.end = None
        self.graph = None
        self.eulergraph = None
        self._endfrag = None
        
        # Creating Double Digests
        dbl1, single_eb1 = RNAGraph._double_digest(self, g_frag, "UC")
        dbl2, single_eb2 = RNAGraph._double_digest(self, uc_frag, "G")
        singles = [e[0] for e in single_eb1 + single_eb2]
        interiors = [e[0] for e in interior_bases(dbl1, dbl2)]

        # Create Graph and Find Graph Start/End Points
        self.graph = to_graph(dbl1, dbl2)
        RNAGraph.find_start(self, singles, interiors)
		
	# Create Eulerian Cycle Version of Graph: Connects Start/End
        self.eulergraph = self.graph.copy()
        try:
            self.eulergraph.remove_edge(self._endfrag[0],self._endfrag[1])
        except (IndexError, nx.exception.NetworkXError):
            pass  # edge may not exist

        end_start_frag = [u for u in self._endfrag] + [self.start]
        edge_label = ''
        if len(end_start_frag) > 2:
            edge_label = ''.join([eb for eb in end_start_frag[1:-1]]) + '*'
        self.eulergraph.add_edge(end_start_frag[0],
                                 end_start_frag[-1],
                                 label=edge_label)
        nodes = list(self.eulergraph.nodes())
        for node in nodes:
            if (not self.eulergraph.out_degree(node)
                and not self.eulergraph.in_degree(node)):
                self.eulergraph.remove_node(node)
        try:
            if self.end != self.start:
                nx.relabel_nodes(self.eulergraph,
                                 {self.end: self.end+'*',
                                  self.start: self.start+'+'},
                                 copy=False)
            else:
                nx.relabel_nodes(self.eulergraph,
                                 {self.end: self.end+'+*'},
                                 copy=False)
        except KeyError:
            pass

    def _double_digest(self, digest_list, split_by):
        """
        Further splits a digest (UC-digest or G-digest) by passed delimiter

        Args:
            digest_list: an iterable digest already split by UC or G enzymes
            split_by: next split ("UC" or "G")

        Returns:
            Tuple (list: extended bases of length 2 or more for graphing,
                   list: single extended bases to determine graph end points)
        """
        new_digest = []
        single_ebs = []
        for base in digest_list:
            if split_by == "UC" and base[-1] != "G":
                self.end = base[-1]
                self._endfrag = to_digest(base, split_by)
            elif split_by == "G" and (base[-1] != "U" and base[-1] != "C"):
                self.end = base[-1]
                self._endfrag = to_digest(base, split_by)

            next_split = to_digest(base, split_by)

            if len(next_split) >= 2:
                new_digest.append(next_split)
            else:
                single_ebs.append(next_split)
        return new_digest, single_ebs

    def find_start(self, singles, interiors):
        """ Find start/end points for RNA sequence """
        diff = list(singles)

        # diff=singles-interiors
        for base in interiors:
            diff.remove(base)

        if diff[0] == self.end:
            self.start = diff[1]
        else:
            self.start = diff[0]
            if not self.end:
                self.end = diff[1]

    def make_seq(self, seq):
        """Returns possible RNA Sequence from Eulerian trail of edges"""
        label = nx.get_edge_attributes(self.graph, 'label')
        final_seq = []
        for (u, v, z) in seq:
            if (label[(u, v, z)] == None):
                if len(final_seq) == 0 or final_seq[-1] != u:
                    final_seq.append(u)
                final_seq.append(v)
            else:
                if len(final_seq) == 0 or final_seq[-1] != u:
                    final_seq.append(u)
                final_seq.append(label[(u, v, z)])
                final_seq.append(v)
        return ''.join(final_seq)

