# sars-covid_google_mat.py
# finds the Google matrix for the PPI with both SARS-human and COVID-human
# interactions

import networkx as nx
import re
import numpy as np
import pandas as pd

# make a networkx graph to hold the human PPI
print('Creating human PPI')
ppi_graph = nx.Graph()

# need the key for converting STRING IDs to uniprot IDs before reading in the
# human PPI
id_dict = dict()
with open('string_to_uniprot_id_dict.tsv', 'r') as in_file:
    # skip header
    next(in_file)
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        # string ID as key, uniprot ID as value
        id_dict[cols[0]] = cols[1]

# now read in the human PPI
with open('string_human_ppi.txt', 'r') as in_file:
    # skip header
    x = next(in_file)
    for line in in_file:
        cols = line.rstrip('\n').split(' ')
        # convert IDs, but since most IDs won't map, use a try/except loop
        try:
            uniprot_ids = [id_dict[cols[0]], id_dict[cols[1]]]
            # filter out low-scoring interactions
            if float(cols[2]) >= 700:
                ppi_graph.add_edge(
                    uniprot_ids[0],
                    uniprot_ids[1]
                )
        except KeyError:
            continue
 
# now add in edges for human-covid interactions
print('Adding human-COVID edges')
# also make a list for later
covid_prots = list()
with open('supp_table_2.tsv', 'r') as in_file:
    # skip both header lines
    x = next(in_file)
    x = next(in_file)
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        # nobody likes whitespace
        covid_prot = re.sub(' ', '_', cols[0])
        human_prot = cols[7]
        covid_prots.append(covid_prot)
        ppi_graph.add_edge(covid_prot, human_prot)

# now add in edges for human-SARS interactions
print('Adding human-SARS edges')
sars_prots = list()
with open('sars_human_interactions.tsv', 'r') as in_file:
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        # uniprot IDs are prefixed by 'uniprotkb:' 
        prot1 = cols[0].lstrip('uniprotkb:')
        prot2 = cols[1].lstrip('uniprotkb:')
        # SARS proteins aren't consistently in the same column
        if prot1.endswith('CVHSA'):
            sars_prots.append(prot1)
        elif prot2.endswith('CVHSA'):
            sars_prots.append(prot2)
        ppi_graph.add_edge(prot1, prot2)

# compute the Google matrix for this network
print('Computing Google matrix')
# make the stochastic matrix from the adjacency matrix
A = nx.adjacency_matrix(ppi_graph)
S = np.zeros(A.shape)
for i in range(A.shape[0]):
    # there shouldn't be any isolated nodes but just in case this is how they
    # would be handled
    if A[i].sum() == 0:
        S[i].fill(1/A.shape[0])
    else:
        # each nonzero entry in this column should be 1/(number of nonzeros)
        # np.nonzero(A[i]) gives a tuple with one array with the indices of the
        # nonzero entries in A[i]
        # since A is an adjacency matrix, the column sum is the same as the
        # number of nonzero elements
        S[i][np.nonzero(A[i])[0]] = 1/A[i].sum()

# make the stochastic matrix into the Google matrix
alpha = 0.85 # I have literally no idea how to choose alpha
G = alpha*S + (1-alpha)*A.shape[0]

print('Saving Google matrix')
np.save('sars-covid_google_mat.npy', G)
