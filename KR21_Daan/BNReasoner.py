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

    def dsep(self,independent1,independent2,given3):
        d_seperated=nx.algorithms.d_separated(self.bn.structure,{independent1},{independent2},{given3})
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

    #Given query variables Q and a possibly empty evidence E, compute the
    # marginal distribution P (Q|E) (12pts). (Note that Q is a subset of the variables in the Bayesian
    # network X with Q ⊂X but can also be Q = X.)
    def marginaldistrubution(self,vars):
        variables = self.variables

        for i in variables:
            cpt= self.bn.get_cpt(i)
            p_positive=0
            p_negative=0

    def deletevariable(self,variable):
        cpts = []
        #get all cpts that mention variable
        for i in self.variables:
            cpt = self.bn.get_cpt(i)
            if variable in cpt.columns:
                cpts.append(cpt)
        print(cpts)

    def grotetabel(self):
        variable = 'Winter?'
        cpt= self.bn.get_cpt(variable)
        for index, row in cpt.iterrows():
            if row['Winter?']== False:
                p_winter_false = row['p']
            if row['Winter?']== True:
                p_winter_true = row['p']

        variable = 'Sprinkler?'
        cpt= self.bn.get_cpt(variable)
        data ={'Sprinkler?': [], 'Winter?': [],'p':[]}
        for index, row in cpt.iterrows():
            if row['Winter?']== False:
                data['p'].append(row['p']*p_winter_false)
                sprinkler = row['Sprinkler?']
                print(sprinkler)
                data['Sprinkler?'].append(row['Sprinkler?'])
                data['Winter?'].append(False)

            elif row['Winter?']== True:
                data['p'].append(row['p'] * p_winter_true)
                data['Sprinkler?'].append(row['Sprinkler?'])
                data['Winter?'].append(True)

        df = pd.DataFrame(data)

        # Print the data
        print(df)


        data = {'Sprinkler?': [], 'Winter?': [],'Rain?':[],'p':[]}
        for index, row in df.iterrows():
            for index2, row2 in self.bn.get_cpt('Rain?').iterrows():
                if row['Winter?'] == row2['Winter?']:
                    data['p'].append(row['p'] * row2['p'])
                    data['Sprinkler?'].append(row['Sprinkler?'])
                    data['Winter?'].append(row['Winter?'])
                    data['Rain?'].append(row2['Rain?'])
        df = pd.DataFrame(data)

        # Print the data
        print(df)




        #multiply all variables that mention variable










bayes = BNReasoner('testing/lecture_example.BIFXML')
#bayes.marginaldistrubution('Wet Grass?')
#bayes.deletevariable('Rain?')
bayes.grotetabel()
# print(bayes.bn.get_all_cpts())
# print(bayes.dsep('Winter?','Slippery Road?','Sprinkler?'))
#print(bayes.MinDegreeOrder())
#print(bayes.MinFillOrder())

