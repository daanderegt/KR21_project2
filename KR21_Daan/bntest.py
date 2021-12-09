import itertools
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
        # print(self.bn.get_all_cpts())
        # G=self.bn.get_interaction_graph()
        # self.bn.draw_structure()
        self.net = net

    # TODO: This is where your methods should go

    def dsep(self, independent1, independent2, given3):
        d_seperated = nx.algorithms.d_separated(self.bn.structure, {independent1}, {independent2}, {given3})
        return d_seperated

    def MinDegreeOrder(self):
        G = self.bn.get_interaction_graph()
        degrees = {}
        for i in self.variables:
            degrees[i] = G.degree(i)
        sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1]))
        return sorted_degrees.keys()

    def MinFillOrder(self):
        degrees = {}
        for i in self.variables:
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

    # Given query variables Q and a possibly empty evidence E, compute the
    # marginal distribution P (Q|E) (12pts). (Note that Q is a subset of the variables in the Bayesian
    # network X with Q ⊂X but can also be Q = X.)
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
        for e in evidence:
            variable = e
            value = evidence[e]
            cpt = self.updating_cpt(variable,value)
            self.bn.update_cpt(variable,cpt)
        print(self.bn.get_all_cpts())
        cpts = []
        for i in variables:
            cpt = self.bn.get_cpt(i)
            cpts.append(self.variable_elimination(cpt))
        cpt1 = cpts[0]
        cpt2 = cpts[1]
        cpt1 = cpt1.merge(cpt2)
        cpt1['p'] = cpt1['p_x'] * cpt1['p_y']
        cpt1 = cpt1.drop(['p_x', 'p_y'], axis=1)
        return cpts

    def variable_elimination(self,cpt):

        for i in cpt.columns[:-2]:
            cpt_next = self.bn.get_cpt(i)
            cpt = self.multiply_cpt(cpt_next, cpt)
            cpt = self.summing_out(cpt,i)
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
        #for i in variables_left:

        # for index, row in cpt.itercolumns():
        #     if row[variables_left] ==
        # for i in cpt[variables_left]:

        cpt = cpt.groupby(variables_left).agg({'p': 'sum'})
        cpt.reset_index(inplace=True)
        return cpt

    def updating_cpt(self,variable,value):
        cpt = self.bn.get_cpt(variable)
        j = 0
        # cpt = cpt.drop('p',axis=1)
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
        return cpt




bayes = BNReasoner('testing/lecture_example.BIFXML')
# print(bayes.bn.get_all_cpts())
# print(bayes.dsep('Winter?','Slippery Road?','Sprinkler?'))
# print(bayes.MinDegreeOrder())
# print(bayes.MinFillOrder())


# print(bayes.marginal_distribution(['Wet Grass?']))
cpt=bayes.marginal_distribution(['Wet Grass?','Slippery Road?'], {'Winter?': True, 'Sprinkler?': True, 'Rain?': True})
print(cpt)
#print(bayes.grotetabel())


# print(bayes.updating_cpt('Wet Grass?', True))