#!/usr/bin/env python3
'''Digest Functions'''
import re
import networkx as nx


def to_digest(rna_seq: object, split_by: object) -> object:
    """Splits an RNA sequence by delimiter enzyme.
    Args:
        rna_seq: an rna sequence string.
        split_by: enzyme delimiter, usually "UC" or "G".
    Returns:
        a list of bases split by provided enzyme delimiter.
    """
    # re: '([^split_by]*[split_by])'
    # '()': keeps delimiter text #[^]:complement set #* :wildcard
    split = re.split('([^' + split_by + ']*[' + split_by + '])', rna_seq)
    return [base for base in split if base]


def input_error(frag1, frag2):
    """ Tests user input error: True if fragmentations do not match """
    f1 = ''.join(frag1)
    f2 = ''.join(frag2)
    return not len(f1) == len(f2)



# Graph Functions

def addv(vertex, graph):
    if vertex not in graph.nodes():
        graph.add_node(vertex)


def dbl_to_graph(dbl, g):
    """ Adds a double-digest (processed by both UC, G enzymes) to NetworkX
    MultiDiGraph
    :param dbl: a double digest (iterable) ex. [['U','G'],['AC','A','G']...]
    :param graph: a networkx multidigraph
    :return: updated networkx multidigraph
    """
    for eb_list in dbl:
        if len(eb_list) == 2:  # for |fragment|==2: ['U']-->['G']
            addv(eb_list[0], g)
            addv(eb_list[1], g)
            g.add_edge(eb_list[0], eb_list[1], label='')
            g.new_edge_key(eb_list[0], eb_list[1])
        else: # for |fragment| > 2: ['U'] -'CAG'-> ['AG'] = node -edgelabel-> node
            addv(eb_list[0], g)
            addv(eb_list[-1], g)
            edge_label = ''.join(eb_list[1:-1])
            g.add_edge(eb_list[0], eb_list[-1], label=edge_label)
            g.new_edge_key(eb_list[0], eb_list[-1])
    return g


def to_graph(dbl1, dbl2):
    """
    Turns 2 double digests of form [['U','C','AG'],['U','G'], ['C','A']]
    into a NetworkX MultiDiGraph (order of parameters doesn't matter)
    :param dbl1: an rna sequence that has been processed by UC enzyme,
                    then G enzyme
    :param dbl2: an rna sequence that has been processed by G enzyme,
                    then UC enzyme
    :return: NetworkX MultiDiGraph
    """
    g = nx.MultiDiGraph()
    dbl_to_graph(dbl1, g)
    dbl_to_graph(dbl2, g)
    return g


# Take 2 double digest lists of form: [['U','C','AG'],['U','G'], ['C','A']]
# Returns a list of interior bases from both digests
def interior_bases(dbl1, dbl2):
    """
    Finds interior bases (middle bases, sandwiched between beginning and end)
    from two double digests
    """
    ib_list = []
    ib1 = [base for base in dbl1 if len(base) > 2]
    ib2 = [base for base in dbl2 if len(base) > 2]
    ib = ib1 + ib2
    for base in ib:
        # interior_list: sublist of form: "[[element],[element]]"
        interior_list = [[e] for e in base[1:-1]]
        ib_list += interior_list
    return ib_list

