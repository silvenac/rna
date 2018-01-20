import re
import pygraphviz as pgv
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import rnagraph as rna
from digest import to_digest
from dfs import custom_edge_dfs


def prep_rnagraph(rna):
    """Creates arguments for RNAGraph class object (for testing)"""
    g_digest = to_digest(rna, "G")
    uc_digest = to_digest(rna, "UC")
    return g_digest, uc_digest


def test_graph(_rna):
    """Tester: Tests graph creation, path finder functions of RNAGraph class
    1. Creates UC, G digests from an RNA seq
    2. Creates graph and euler graph from digests
    3. Uses graph to find Eulerian Trails
    4. Uses euler graph to save png of Eulerian graph representation of RNA seq
    """
    g,uc = prep_rnagraph(_rna)
    rg = rna.RNAGraph(g,uc)
    edges = list(rg.graph.edges(keys=True))
    nodes = list(rg.graph.nodes())

    allpaths = custom_edge_dfs(rg.graph,rg.start)
    joinedpaths = set([rg.make_seq(path) for path in allpaths])
    print("Original Sequence:",_rna)
    print("G:",g)
    print("UC:",uc)
    A = nx.nx_agraph.to_agraph(rg.eulergraph)
    A.layout(prog='dot')
    A.draw(_rna+'.png')
    print("Start:",rg.start,"End:",rg.end)
    print("Paths")
    for i,path in enumerate(joinedpaths):
        print(i,path)
    print()


def main():
    test_graph("AUGAUCGGACUAUACGCU")
    test_graph("AGUCAGUGAGCA")
    test_graph("GCAGAAAAAAACCUUAAGUCUGC")
    test_graph("GGCUGUUACCGAAAAAGCAGCAGGCAGAGC")


if __name__ == '__main__':
    main()
