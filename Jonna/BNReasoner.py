from copy import deepcopy
import itertools
from typing import Union

import networkx

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

    def d_sep_networkx(self,independent1,independent2,given3): # To verify if our own method gives the same answer
        d_seperated=nx.algorithms.d_separated(self.bn.structure,{independent1},{independent2},{given3})
        return d_seperated

    def MinDegreeOrder(self, variables, order):
        G = self.bn.structure.to_undirected()
        degrees = {}
        all_combinations = []
        neighbors = []
        for variable in variables:
            degrees[variable] = len(list(nx.neighbors(G, variable)))

        least_edges = min(degrees, key=degrees.get)
        for i in nx.neighbors(G, least_edges):
            neighbors.append(i)

        for i in range(len(list(neighbors)) + 1):
            neighbors_iterated = itertools.combinations(neighbors, i)
            combinations_list = list(neighbors_iterated)
            all_combinations += combinations_list

        neighbors = [i for i in all_combinations if len(i) == 2]
        for neighbor in neighbors:
            self.bn.add_edge(neighbor)
            print(neighbor)
        order.append(least_edges)
        self.bn.del_var(least_edges)
        variables = self.bn.get_all_variables()
        self.bn.draw_structure()

        if len(variables) == 0:
            return order
        else:
            return self.MinDegreeOrder(variables, order)

    def MinFillOrder(self):
        degrees = {}
        for i in self.variables:
            degrees[i] = self.check_edges_del_var(i)
        sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1], reverse=True))
        return sorted_degrees.keys()

    def MinFillOrderr(self, variables, order):
        G = self.bn.structure.to_undirected()
        fill = {}
        for variable in variables:
            length=len(self.getneighborsedges(variable))
            fill[variable]=length

        sorted_degrees = dict(sorted(fill.items(), key=lambda item: item[1]))
        max=list(sorted_degrees.keys())[0]
        order.append(max)
        self.del_variable(G,max)
        variables = self.bn.get_all_variables()
        if len(variables) == 0:
            return order
        else:
            return self.MinFillOrderr(variables, order)



    def getneighborsedges(self,node):
        G = self.bn.structure.to_undirected()

        neighbors=[]
        all_combinations=[]
        for i in nx.neighbors(G, node):
            neighbors.append(i)

        for i in range(len(list(neighbors)) + 1):
            neighbors_iterated = itertools.combinations(neighbors, i)
            combinations_list = list(neighbors_iterated)
            all_combinations += combinations_list
        neighbors = [i for i in all_combinations if len(i) == 2]
        return neighbors

    def del_variable(self,G,variable):
        edges=self.getneighborsedges(variable)
        for edge in edges:
            self.bn.add_edge(edge)
        self.bn.del_var(variable)

    def check_edges_del_var(self, var):
        bay = BayesNet()
        # Loads the BN from an BIFXML file
        bay.load_from_bifxml(self.net)
        bay.del_var(var)
        edges = bay.get_number_of_edges()
        return edges

    def node_pruning(self, Q_variables, E_evidence):
        G = self.bn.get_digraph()
        for i in self.variables:
            if i not in Q_variables and i not in E_evidence and G.out_degree(i) == 0:
                self.bn.del_var(i)

        self.bn.draw_structure()

        if E_evidence:
            for e in E_evidence:
                evidence = e
                booleanvalue = E_evidence[e]
                self.edge_pruning(evidence, booleanvalue)

    def edge_pruning(self, evidence, booleanvalue):
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

    ''' Marginal distribution '''
    def multiply_cpt(self,cpt1, cpt2):
        cpt1_list=list(cpt1.columns)
        cpt2_list = list(cpt2.columns)
        common_vars = [var for var in cpt1_list if var in cpt2_list and var != 'p']
        if not common_vars:
            return 'Multiplication is not possible because there are no common variables in the given CPTs'

        cpt1 = cpt1.merge(cpt2,left_on=common_vars,right_on=common_vars)
        cpt1['p']= cpt1['p_x']*cpt1['p_y']
        cpt1 = cpt1.drop(['p_x','p_y'],axis=1)
        return cpt1

    def summing_out(self, cpt , variable):
        #als andere variabelen gelijk zijn in cpt moeten deze vermenigvuldigd worden
        cpt = cpt.drop([variable],axis=1)
        variables_left = [var for var in cpt.columns if var != variable and var != 'p']
        cpt = cpt.groupby(variables_left).agg({'p': 'sum'})
        cpt.reset_index(inplace=True)
        return cpt

    def variable_elimination(self,cpt):
        order = self.MinFillOrder(cpt.columns[:-2])
        length = len(cpt.columns)
        if length == 2:
            return cpt
        order1 = cpt.columns[:-2]
        for i in order1:
            cpt_next = self.bn.get_cpt(i)
            cpt = self.multiply_cpt(cpt_next, cpt)
            print(cpt)
            cpt = self.summing_out(cpt,i)

        return self.variable_elimination(cpt)

    def updating_cpt(self,variable,value):
        cpt = self.bn.get_cpt(variable)
        j = 0
        for i in cpt[variable]:
            if value:
                if i:
                    cpt['p'][j] = 1
                else:
                    cpt['p'][j] = 0
            else:
                if i:
                    cpt['p'][j] = 0
                else:
                    cpt['p'][j] = 1
            j += 1
        print('cpt=',variable,cpt)
        return cpt

    def marginal_distribution(self,variables,evidence):
        if evidence:
            for e in evidence:
                variable = e
                value = evidence[e]
                cpt = self.updating_cpt(variable,value)
                self.bn.update_cpt(variable,cpt)
        cpts = []
        for i in variables:
            cpt = self.bn.get_cpt(i)
            cpts.append(self.variable_elimination(cpt))

        cpt1 = cpts[0]
        for i in cpts[1:]:
            Data = {}
            cpt2 = i
            for i in cpt1:
                if i != 'p':
                    Data.update({i: []})
            for j in cpt2:
                if j != 'p':
                    Data.update({j: []})
            Data.update({'p':[]})
            for index, row in cpt1.iterrows():
                for index2, row2 in cpt2.iterrows():
                    for i in Data.keys():
                        if i != 'p':
                            if i in cpt1:
                                Data[i].append(row[i])
                            if i in cpt2:
                                Data[i].append(row2[i])
                        else:
                            Data['p'].append(row['p'] * row2['p'])
            cpt1 = pd.DataFrame(Data)
        return cpt1


    ''' MAP '''

    def MAP(self, variables, evidence):
        cpt = self.marginal_distribution(variables,evidence)
        max = cpt['p'].max()
        maxrow=cpt.loc[cpt['p']==max]
        return maxrow

    ''' MPE '''

    def MPE(self,variables,evidence):
        evidence = {}
        maxrow = self.MAP(variables,evidence)
        return maxrow


# !! self.bn.get_interaction_graph() = self.bn.get_digraph()

if __name__ == "__main__":
    bayes = BNReasoner('testing/corona_example.BIFXML')
    #bayes = BNReasoner('testing/b500-31.xml')
    #bayes.node_pruning('node100', [('node203', True), ('node333', False), ('node1', False), ('node33', False), ('node400', False)])
    #bayes.d_sep('node5', 'node30', ['node22', 'node2'])
    #print(bayes.marginal_distribution(['Corona?'], {'Obesitas?': True}))
    #print (bayes.marginal_distribution(['node22'],{'node333': False}))
    #print(bayes.node_pruning(['node22'],{'node333': False}))

    # print(bayes.MinDegreeOrder(bayes.bn.get_all_variables(), []))
    #bayes.del_variable(bayes.bn.structure.to_undirected(),'Winter?')
    print(bayes.MinFillOrderr(bayes.variables,[]))
    # print(bayes.MinFillOrder())

