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
            self.bn = BayesNet()
            self.bn.load_from_bifxml(net)

        else:
            self.bn = net
        self.variables= self.bn.get_all_variables()
        self.net = net

    def d_sep(self,var1,var2,given3): 
        
        graph = self.bn.get_interaction_graph()
        for node in given3: 
            graph.remove_node(node)
        
        try:
            sub = graph.subgraph(nx.shortest_path(graph.to_undirected(), var1, var2)) 
            print(False)
            return False
        except:
            print(True)
            return True

    def dsep(self,independent1,independent2,given3): # To verify if our own method gives the same answer
        d_seperated=nx.algorithms.d_separated(self.bn.structure,{independent1},{independent2},{given3})
        return d_seperated

    ''' Add MinDegreeOrder '''

    ''' Add MinFillOrder '''

    ''' Add check_edges_del_var '''

    def node_pruning(self, Q_variables, E_evidence):
        G = self.bn.get_digraph()
        for i in self.variables:
            if i not in Q_variables and i not in E_evidence and G.out_degree(i) == 0:
                self.bn.del_var(i)

        self.bn.draw_structure()
        for evidence in E_evidence:
            self.edge_pruning(evidence)

    def edge_pruning(self, given_evidence):
        evidence = given_evidence[0]
        booleanvalue = given_evidence[1]

        G = self.bn.get_digraph()

        for i in self.bn.get_all_variables():
            if i in evidence and G.out_degree(i) > 0:
                edges = G.out_edges(i)
                for y in edges:
                    try: 
                        self.bn.del_edge(y)
                    except:
                        continue

        self.bn.draw_structure()

        cpts_evidence = []
        for i in self.bn.get_all_variables():
            cpt = self.bn.get_cpt(i)

            if evidence in cpt.columns[0:-2]:
                cpts_evidence.append(cpt)

        for i in cpts_evidence:
            cpt_evidence = self.bn.get_compatible_instantiations_table(pd.Series({evidence: booleanvalue}), i)

            df = pd.DataFrame(cpt_evidence)
            df.drop(df[df[evidence] != booleanvalue].index, inplace = True)

            del df[evidence]
            print("Modified CPT's: ", df)

    ''' Marginal distributions '''

    ''' MAP '''

    ''' MPE '''


# !! self.bn.get_interaction_graph() = self.bn.get_digraph()

if __name__ == "__main__":
    bayes = BNReasoner('Jonna/testing/b500-31.xml')
    #bayes.node_pruning('node100', [('node203', True), ('node333', False), ('node1', False), ('node33', False), ('node400', False)])
    #bayes.d_sep('node5', 'node30', ['node22', 'node2'])

