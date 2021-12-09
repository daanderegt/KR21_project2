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

    def marginal_distribution(self,variables):
        for variable in variables:
            cpt = self.bn.get_cpt(variable)
            for i in cpt.columns[:-2]:
                if not i:
                    return 'Tabel'
                else:
                    cpt_next = self.bn.get_cpt(i)
                    cpt = self.multiply_cpt(cpt, cpt_next)



    def J_P_D(self, variables):
        cpt = bayes.bn.get_cpt(variables[1])
        for i in variables[1:]:
            cpt = bayes.multiply_cpt(cpt, bayes.bn.get_cpt(i))
        return cpt

    def summing_out(self, cpt , variable):

        return 'iets'


    def variable_elimination(self,variable):
        #get all cpts that mention variable
        variables=[]
        for i in self.variables:
            cpt=self.bn.get_cpt(i)
            if variable in cpt.columns:
                variables.append(cpt.columns[-2])

        #get Joint probability distrubution of these variables
        jpd = self.J_P_D(variables)
        self.summing_out(jpd,variable)








bayes = BNReasoner('testing/lecture_example.BIFXML')
# print(bayes.bn.get_all_cpts())
# print(bayes.dsep('Winter?','Slippery Road?','Sprinkler?'))
# print(bayes.MinDegreeOrder())
# print(bayes.MinFillOrder())


print(bayes.marginal_distribution(['Wet Grass?']))
#print(bayes.grotetabel())
