from typing import Union
from BayesNet import BayesNet
import networkx as nx
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
        self.bn.draw_structure()
        self.net = net

    # TODO: This is where your methods should go

    def dsep(self,independent1,independent2,given3):
        G = self.bn.get_interaction_graph()
        d_seperated=nx.d_separated(G,independent1,independent2,given3)
        return d_seperated

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



bn = BNReasoner('testing/lecture_example.BIFXML')
#print(bn.dsep('Winter?','Slippery Road?','Sprinkler?'))
print(bn.MinDegreeOrder())
print(bn.MinFillOrder())
