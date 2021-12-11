import datetime
from typing import Union
from BayesNet import BayesNet
import networkx as nx
import pandas as pd
import numpy as np
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
        self.variables = self.bn.get_all_variables()
        self.net = net

    # TODO: This is where your methods should go

    def dsep(self, independent1, independent2, given3):
        d_seperated = nx.algorithms.d_separated(self.bn.structure, {independent1}, {independent2}, {given3})
        return d_seperated

    def MinDegreeOrder(self,variables):
        G = self.bn.get_interaction_graph()
        degrees = {}
        for i in variables:
            degrees[i] = G.degree(i)
        sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1]))
        return sorted_degrees.keys()

    def MinFillOrder(self, variables):
        degrees = {}
        for i in variables:
            degrees[i] = self.check_edges_del_var(i)
        sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1], reverse=True))
        return sorted_degrees.keys()

    def check_edges_del_var(self, var):
        bay = BayesNet()
        # Loads the BN from an BIFXML file
        bay.load_from_bifxml(self.net)
        bay.del_var(var)
        edges = bay.get_number_of_edges()
        return edges

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

    def variable_elimination(self,cpt):
        order = self.MinFillOrder(cpt.columns[:-2])
        print(order)
        order1 = cpt.columns[:-2]
        for i in order:
            cpt_next = self.bn.get_cpt(i)
            cpt = self.multiply_cpt(cpt_next, cpt)
            print(cpt)
            cpt = self.summing_out(cpt,i)
            print(i)
            print(cpt)
            if cpt_next.columns[:-2].empty:
                return cpt
        return self.variable_elimination(cpt)


    def J_P_D(self, variables):
        cpt = bayes.bn.get_cpt(variables[1])
        for i in variables[1:]:
            cpt = bayes.multiply_cpt(cpt, bayes.bn.get_cpt(i))
        return cpt

    def summing_out(self, cpt , variable):
        #als andere variabellen gelijk zijn in cpt moeten deze vermenigvuldigd worden
        cpt = cpt.drop([variable],axis=1)
        variables_left = [var for var in cpt.columns if var != variable and var != 'p']
        cpt = cpt.groupby(variables_left).agg({'p': 'sum'})
        cpt.reset_index(inplace=True)
        return cpt

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

    def MAP(self, variables, evidence):
        cpt = self.marginal_distribution(variables,evidence)
        max = cpt['p'].max()
        maxrow=cpt.loc[cpt['p']==max]
        return maxrow

    def MPE(self,variables,evidence):
        evidence = {}
        maxrow = self.MAP(variables,evidence)
        return maxrow



bayes = BNReasoner('testing/lecture_example.BIFXML')
# print(bayes.bn.get_all_cpts())
# print(bayes.dsep('Winter?','Slippery Road?','Sprinkler?'))
# print(bayes.MinDegreeOrder())
# print(bayes.MinFillOrder())
# print(bayes.marginal_distribution(['Wet Grass?']))
# cpt=bayes.marginal_distribution(['Wet Grass?','Slippery Road?','Rain?'], {'Winter?': True, 'Sprinkler?': False })
# print(cpt)
#print(bayes.MAP(['Wet Grass?','Slippery Road?','Rain?'], {'Winter?': True, 'Sprinkler?': False }))execution_time = timeit.timeit(code, number=1)
begin_time = datetime.datetime.now()
print(bayes.MPE(['Wet Grass?'], {'Winter?': True, 'Sprinkler?': False }))
print(datetime.datetime.now() - begin_time)
