from copy import deepcopy
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
        
        graph = self.bn.get_interaction_graph() # create interaction graph
        for node in given3: # remove evidence nodes and edges
            graph.remove_node(node)

        # nx.draw(graph, with_labels = True)
        # plt.show()
        
        try:
            sub = graph.subgraph(nx.shortest_path(graph.to_undirected(), var1, var2)) # check if there still exists a path between var1 and var2
            return False
        except:
            return True

    def dsep(self,independent1,independent2,given3): # To verify if our own method gives the same answer
        d_seperated=nx.algorithms.d_separated(self.bn.structure,{independent1},{independent2},{given3})
        return d_seperated

    # def MAP(self, dict): #placeholder
    #     df = compute_posterior_marginal(dict)
    #     for column in df:
    #         answer = df["p"].idmax()

    #     return answer

    # def MPE(self, table):
    #     answer = table["p"].idmax()


    #     return answer

    def multiply(self, factor1, factor2):
        var1 = factor1.columns[-2]
        var2 = factor2.columns[-2]

        F2F1 = deepcopy(factor2)

        for i in range(0, len(factor1[var1])):
            if factor1[var1].values[i] == True:
                True_p = factor1['p'][i]
            if factor1[var1].values[i] == False:
                False_p = factor1['p'][i]

        for i in range(0, len(factor2[var1])):
            if factor2[var1].values[i] == True:
                F2F1.at[i, 'p'] = factor2['p'][i] * False_p
            if factor2[var1].values[i] == False:
                F2F1.at[i, 'p'] = factor2['p'][i] * False_p

        return F2F1

       
    def summing_out(self, params):
        table = "table"

        return table 

bayes = BNReasoner('Jonna/testing/lecture_example.BIFXML')

# print(bayes.dsep('Slippery Road?', 'Wet Grass?', 'Rain?'))
# print(bayes.d_sep('Slippery Road?', 'Winter?', ['Rain?', 'Wet Grass?']))

cpts = bayes.bn.get_all_cpts()
print(bayes.multiply(cpts["Winter?"], cpts["Sprinkler?"]))


