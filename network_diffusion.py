# network_diffusion.py
# take the COVID-19-human interactions, join them to a "normal" human PPI, and
# do network diffusion

import networkx as nx
import re
import scipy as sp
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
print('Adding human-virus edges')
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

# now add in edges for human-SARS interactinos
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

# compute diffusion matrix
print('Computing Laplacian matrix')
lap_mat = nx.laplacian_matrix(ppi_graph).tocsc()
id_mat = sp.sparse.identity(lap_mat.shape[0]).tocsc()
print('Computing diffusion matrix')
diff_mat = sp.sparse.linalg.inv(id_mat + (0.1 * lap_mat))

print('Doing COVID-other diffusion')
# need to create a binary vector with 1s for all COVID proteins and 0s for all
# human genes in the same order as the rows in diff_mat
covid_diff_vec = list()
# testing on small graphs indicates that the order of nodes in graph.nodes()
# always matches the order of columns in the Laplacian matrix
for prot in ppi_graph.nodes():
    if prot in covid_prots:
        covid_diff_vec.append(1)
    else:
        covid_diff_vec.append(0)

covid_diff_vec = np.array(covid_diff_vec)

# now multiply the diffusion matrix with this vector to get our ranked list of
# proteins
covid_diff_result = diff_mat * covid_diff_vec

print('Doing SARS-other diffusion')
# need to create a binary vector with 1s for all COVID genes and 0s for all
# human genes in the same order as the rows in diff_mat
sars_diff_vec = list()
# testing on small graphs indicates that the order of nodes in graph.nodes()
# always matches the order of columns in the Laplacian matrix
for prot in ppi_graph.nodes():
    if prot in sars_prots:
        sars_diff_vec.append(1)
    else:
        sars_diff_vec.append(0)

sars_diff_vec = np.array(sars_diff_vec)

# now multiply the diffusion matrix with this vector to get our ranked list of
# proteins
sars_diff_result = diff_mat * sars_diff_vec

print('Preparing output')
# doesn't come with labels but we know the order of the elements in diff_result
# and ppi_graph.nodes() is the same
covid_df = pd.DataFrame(
    list(zip(ppi_graph.nodes(), covid_diff_result)),
    columns = ['protein', 'rank']
)
covid_df.to_csv('covid_diff.csv', index = False)

# doesn't come with labels but we know the order of the elements in diff_result
# and ppi_graph.nodes() is the same
sars_df = pd.DataFrame(
    list(zip(ppi_graph.nodes(), sars_diff_result)),
    columns = ['protein', 'rank']
)
sars_df.to_csv('sars_diff.csv', index = False)
