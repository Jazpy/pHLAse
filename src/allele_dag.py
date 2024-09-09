'''
Directed Acyclic Graph representing possible haplotype phases
'''
import math

import networkx as nx
import matplotlib.pyplot as plt


class ADAG:
    def __init__(self):
        self.dag = nx.DiGraph()
        self.dag.add_node(0, allele='x')
        self.last_level = [0]
        self.nodes      = 1


    def add_level(self, linkages):
        fsts = [l[0][0] for l in linkages]
        snds = [l[0][1] for l in linkages]
        lnks = [l[1] for l in linkages]

        # Get log likelihoods
        wgts = [math.log(l / sum(lnks)) for l in lnks]

        preds = []
        for n in self.last_level:
            preds.extend(self.dag.predecessors(n))

        for f in fsts:
            found = False
            for n in self.last_level:
                if self.dag.nodes(data=True)[n]['allele'] == f:
                    found = True
                    break
            if not found:
                self.dag.add_node(self.nodes, allele=f)
                for p in preds:
                    self.dag.add_edge(p, self.nodes, weight=1, links=0)
                self.last_level.append(self.nodes)
                self.nodes += 1

        new_level = []
        allele_to_id = {}
        for s in set(snds):
            self.dag.add_node(self.nodes, allele=s)
            allele_to_id[s] = self.nodes
            new_level.append(self.nodes)
            self.nodes += 1

        for n in self.last_level:
            added = False
            for f, s, l, w in zip(fsts, snds, lnks, wgts):
                if self.dag.nodes(data=True)[n]['allele'] == f:
                    self.dag.add_edge(n, allele_to_id[s], weight=w, links=l)
                    added = True

            if not added:
                for s in snds:
                    self.dag.add_edge(n, allele_to_id[s], weight=1, links=0)

        self.last_level = new_level


    def add_unlinked_level(self, alleles):
        new_level = []
        for n in alleles:
            self.dag.add_node(self.nodes, allele=n)

            for m in self.last_level:
                self.dag.add_edge(m, self.nodes, weight=1, links=0)

            new_level.append(self.nodes)
            self.nodes += 1

        self.last_level = new_level


    def get_paths(self):
        ret = []
        for path in nx.all_simple_paths(self.dag, 0, self.last_level):
            ret.append([self.dag.nodes(data=True)[x]['allele'] for x in list(path)[1:]])
        return ret


    def draw(self):
        pos = nx.bfs_layout(self.dag, start=0)
        labels = nx.get_node_attributes(self.dag.subgraph(list(self.dag.nodes)[1:]), 'allele')
        weights = nx.get_edge_attributes(self.dag.subgraph(list(self.dag.nodes)[1:]), 'links')
        nx.draw(self.dag.subgraph(list(self.dag.nodes)[1:]), pos, labels=labels, with_labels=True)
        nx.draw_networkx_edge_labels(self.dag, pos, edge_labels=weights)
        plt.show()

