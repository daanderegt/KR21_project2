import itertools
from typing import Union
from BayesNet import BayesNet
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt


class BNReasoner:
    def __init__(self, net: Union[str, BayesNet]):
        """
        :param net: either file path of the bayesian network in BIFXML format or BayesNet object
        """
        if type(net) == str:
            # constructs a BN object
            self.bn = BayesNet()
            # Loads the BN from an BIFXML file
            self.bn.load_from_bifxml(net)

        else:
            self.bn = net
        self.variables= self.bn.get_all_variables()
        #print(self.bn.get_all_cpts())
        #G=self.bn.get_interaction_graph()
        #self.bn.draw_structure()
        self.net = net

    # TODO: This is where your methods should go


    def MinDegreeOrder(self):
        G = self.bn.get_interaction_graph()
        degrees = {}
        for i in self.variables:
            degrees[i]=G.degree(i)
        sorted_degrees=dict(sorted(degrees.items(), key=lambda item: item[1]))
        return sorted_degrees.keys()

    def MinFillOrder(self):
        degrees = {}
        print()
        for i in self.variables:
            degrees[i]=self.check_edges_del_var(i)
        sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1],reverse=True))
        return sorted_degrees.keys()

    def check_edges_del_var(self,var):
        bay = BayesNet()
        # Loads the BN from an BIFXML file
        bay.load_from_bifxml(self.net)
        bay.del_var(var)
        edges = bay.get_number_of_edges()
        return edges

    def is_closed(self, var1, var2):
        return ''

    def check_neighbors(self, given3):
        neighbors = []
        
        graph = self.bn.get_interaction_graph()
        nx.draw(graph, with_labels = True)
        plt.show()

        for neighbor in graph.neighbors(given3):
            neighbors.append(neighbor)
        return neighbors


    def d_sep(self,var1,var2,given3):
        
        graph = self.bn.get_interaction_graph()
        graph.remove_node(given3)

        # nx.draw(graph, with_labels = True)
        # plt.show()
        
        try:
            sub = graph.subgraph(nx.shortest_path(graph.to_undirected(), var1, var2))
            return False
        except:
            return True

    def dsep(self,independent1,independent2,given3):
        d_seperated=nx.algorithms.d_separated(self.bn.structure,{independent1},{independent2},{given3})
        return d_seperated


bayes = BNReasoner('Jonna/testing/lecture_example.BIFXML')

print(bayes.dsep('Slippery Road?', 'Wet Grass?', 'Rain?'))
print(bayes.d_sep('Slippery Road?', 'Wet Grass?', 'Rain?'))

#print(bayes.check_neighbors('Rain?'))
