from typing import Union

import pandas as pd

from BayesNet import BayesNet
import networkx as nx
import matplotlib.pyplot


class BNReasonerChantal:
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

    def MinDegreeOrder(self, variables):
        G = self.bn.get_interaction_graph()
        degrees = {}
        for i in variables:
            degrees[i]=G.degree(i)
        sorted_degrees=dict(sorted(degrees.items(), key=lambda item: item[1]))
        return sorted_degrees

    def MinFillOrder(self, variables):
        degrees = {}
        print()
        for i in variables:
            degrees[i]=self.check_edges_del_var(i)
        sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1],reverse=True))
        return sorted_degrees

    def check_edges_del_var(self,var):
        bay = BayesNet()
        # Loads the BN from an BIFXML file
        bay.load_from_bifxml(self.net)
        bay.del_var(var)
        edges = bay.get_number_of_edges()
        return edges

    ''' Remove leaf nodes that are not part of the given variables and evidence '''
    def node_pruning(self, Q_variables, E_evidence):
        ''' Load interaction graph '''
        G = self.bn.get_interaction_graph()
        ''' For all variables in the interaction graph '''
        for i in self.variables:
            ''' If the variable is not in the given variables and evidence, and the variable has out_degree 0, 
                then delete that variable from the graph'''
            if i not in Q_variables and i not in E_evidence and G.out_degree(i) == 0:
                self.bn.del_var(i)

        ''' Show the new interaction graph with deleted nodes'''
        self.bn.draw_structure()
        for evidence in E_evidence:
            self.edge_pruning(evidence)

    ''' Remove edges from the evidence node, and replace CPT for that evidence node with a smaller CPT '''
    def edge_pruning(self, given_evidence):
        evidence = given_evidence[0]
        booleanvalue = given_evidence[1]

        ''' Load interaction graph '''
        G = self.bn.get_interaction_graph()

        ''' For all variables in the interaction graph '''
        for i in self.bn.get_all_variables():
            ''' If the variable is part of variable in the evidence, and the variable has out_degree > 0, 
                then delete edges from that variable '''
            if i in evidence and G.out_degree(i) > 0:
                edges = G.out_edges(i)
                for y in edges:
                    self.bn.del_edge(y)

        ''' Show the new interaction graph with deleted edges'''
        self.bn.draw_structure()

        ''' Replace CPT for the evidence nodes with smaller CPT: '''

        ''' Get all CPT's where evidence is in the cpt (but not in the last two columns) '''
        cpts_evidence = []
        for i in self.bn.get_all_variables():
            cpt = self.bn.get_cpt(i)

            if evidence in cpt.columns[0:-2]:
                cpts_evidence.append(cpt)

        ''' For each CPT, delete rows and evidence column '''
        for i in cpts_evidence:
            cpt_evidence = self.bn.get_compatible_instantiations_table(pd.Series({evidence: booleanvalue}), i)
            #print("CPT with evidence = ", booleanvalue, ": ", cpt_evidence)

            ''' Delete the evidence rows where evidence is not the booleanvalue (true or false)'''
            df = pd.DataFrame(cpt_evidence)
            df.drop(df[df[evidence] != booleanvalue].index, inplace = True)
            #print("Smaller CPT with rows deleted: ", df)

            ''' Delete the evidence column to end up with the smaller CPT '''
            del df[evidence]
            #print("Smaller CPT with evidence column deleted: ", df)
            print("Modified CPT's: ", df)

bn = BNReasonerChantal('testing/lecture_example.BIFXML')

# print(bn.dsep('Winter?','Slippery Road?','Sprinkler?'))
print(bn.MinDegreeOrder())
print(bn.MinFillOrder())

print(bn.node_pruning('Wet Grass?', [('Winter?', True), ('Rain?', False)]))