""""
BNReasoner
authors: Pleuntje Binnerts, Jonna van Doorn, Daan de Regt, Chantal Vogels
"""

import itertools
from typing import Union
import time
import numpy as np
from BayesNet import BayesNet
import networkx as nx
import pandas as pd

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

    def d_sep(self, var1, var2, given3):
        graph = self.bn.get_interaction_graph()
        for node in given3:
            graph.remove_node(node)

        try:
            for node1 in var1:
                for node2 in var2:
                    graph.subgraph(nx.shortest_path(graph.to_undirected(), node1, node2))
            return False
        except:
            return True

    def MinDegreeOrder(self, variables, order,bn):
        G = bn.structure.to_undirected()
        degrees = {}
        for variable in variables:
            degrees[variable] = len(list(nx.neighbors(G, variable)))

        least_edges = min(degrees, key=degrees.get)
        neighbors = self.getneighborsedges(least_edges,bn)
        for neighbor in neighbors:
            bn.add_edge(neighbor)
        order.append(least_edges)
        bn.del_var(least_edges)

        variables = list(variables)
        variables.remove(least_edges)

        if len(variables) == 0:
            return order
        else:
            return self.MinDegreeOrder(variables, order,bn)

    def MinFillOrder(self, variables, order,bn):
        bn.structure.to_undirected()
        fill = {}
        for variable in variables:
            length=len(self.getneighborsedges(variable,bn))
            fill[variable]=length

        sorted_degrees = dict(sorted(fill.items(), key=lambda item: item[1]))
        max=list(sorted_degrees.keys())[0]
        order.append(max)
        self.del_variable(max,bn)
        variables = list(variables)
        variables.remove(max)
        if len(variables) == 0:
            return order
        else:
            return self.MinFillOrder(variables, order,bn)

    def getneighborsedges(self,node,bn):
        G = bn.structure.to_undirected()
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

    def del_variable(self,variable,bn):
        edges=self.getneighborsedges(variable,bn)
        for edge in edges:
            bn.add_edge(edge)
        bn.del_var(variable)

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
        cpt = cpt.drop([variable],axis=1)
        variables_left = [var for var in cpt.columns if var != variable and var != 'p']
        cpt = cpt.groupby(variables_left).agg({'p': 'sum'})
        cpt.reset_index(inplace=True)
        return cpt

    def variable_elimination(self,cpt,order,function):
        length = len(cpt.columns)
        newbn = BayesNet()
        newbn.load_from_bifxml(self.net)
        if length == 2 or length == len(self.variables)+1:
            return cpt
        if order == 'standard':
            order1 = cpt.columns[:-2]
        elif order == 'degree':
            order1 = self.MinDegreeOrder(cpt.columns[:-2], [],newbn)
        elif order == 'fill':
            order1 = self.MinFillOrder(cpt.columns[:-2], [],newbn)

        for i in order1:
            cpt_next = self.bn.get_cpt(i)
            cpt = self.multiply_cpt(cpt_next, cpt)
            if not function =='MPE':
                cpt = self.summing_out(cpt,i)

        return self.variable_elimination(cpt, order, function)

    def updating_cpt(self,variable,value):
        cpt = self.bn.get_cpt(variable)
        j = 0
        for i in cpt[variable]:
            if value:
                if i == False:
                    cpt = cpt.drop(j)
            else:
                if i:
                    cpt = cpt.drop(j)
            j += 1
        cpt = cpt.reset_index(drop=True)
        return cpt

    def evidence_normalisation(self, variable, value):
        function = ''
        cpt = self.marginal_distribution([variable],{}, function)
        j = 0
        for i in cpt[variable]:
            if i == value:
                p = cpt.loc[j, 'p']
                return p
            j+=1

    def normalisation_factor(self,evidence):
        p = []
        norm_factor = 1
        if evidence:
            for e in evidence:
                variable = e
                value = evidence[e]
                cpt_evidence = (self.evidence_normalisation(variable, value))
                p.append(cpt_evidence)
            for i in p:
                norm_factor = norm_factor * i
        return norm_factor

    def marginal_distribution(self,variables,evidence,order, function):
        cpts = []
        if evidence:
            for e in evidence:
                variable = e
                value = evidence[e]
                cpt = self.updating_cpt(variable,value)
                self.bn.update_cpt(variable,cpt)
        for i in variables:
            cpt = self.bn.get_cpt(i)
            cpts.append(self.variable_elimination(cpt,order, function))

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

    def mainfunction(self,variables,evidence,order, function):
        cpt = self.marginal_distribution(variables,evidence,order, function)
        normalisation_factor = self.normalisation_by_one(cpt)
        for i in cpt['p']:
            cpt = cpt.replace(i, i/normalisation_factor)
        return cpt

    def normalisation_by_one(self,cpt):
        norm_factor = 0
        for i in cpt['p']:
            norm_factor = norm_factor + i
        return norm_factor

    def MAP(self, variables, evidence, order):
        function = 'MAP'
        cpt = self.mainfunction(variables,evidence,order, function)
        max = cpt['p'].max()
        maxrow=cpt.loc[cpt['p']==max]
        return maxrow

    def MPE(self, evidence,order):
        variables = self.variables
        MPE = self.MAP(variables, evidence, order)
        return MPE

if __name__ == "__main__":
    bayes = BNReasoner('testing/lecture_example.BIFXML')
    total = []
    for i in range(1):
        start = time.time()
        print(bayes.MPE({'Winter?': True, 'Symptoms?': True, 'Senior?': True}, 'degree'))
        print(bayes.MAP(['Corona?'], {'Winter?': True, 'Symptoms?': True, 'Senior?': True}, 'standard'))
        end = time.time()
        count_time = round(((end - start) * 1000),4)
        total.append(count_time)

    print(total)
    print('avg:',np.mean(total))
    print('std:',np.std(total))