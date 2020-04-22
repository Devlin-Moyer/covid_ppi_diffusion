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
print('Adding human-COVID edges')
# also make a list for later
interactors = list()
with open('supp_table_2.tsv', 'r') as in_file:
    # skip both header lines
    x = next(in_file)
    x = next(in_file)
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        # nobody likes whitespace
        covid_prot = re.sub(' ', '_', cols[0])
        human_prot = cols[7]
        interactors.append(human_prot)
        ppi_graph.add_edge(covid_prot, human_prot)

# compute diffusion matrix
print('Computing Laplacian matrix')
lap_mat = nx.laplacian_matrix(ppi_graph).tocsc()
id_mat = sp.sparse.identity(lap_mat.shape[0]).tocsc()
print('Computing diffusion matrix')
diff_mat = sp.sparse.linalg.inv(id_mat + (0.1 * lap_mat))

print('Creating diffusion vector')
# need to create a binary vector with 1s for all COVID genes and 0s for all
# human genes in the same order as the rows in diff_mat
diff_vec = list()
# testing on small graphs indicates that the order of nodes in graph.nodes()
# always matches the order of columns in the Laplacian matrix
for prot in ppi_graph.nodes():
    if prot in interactors:
        diff_vec.append(1)
    else:
        diff_vec.append(0)

diff_vec = np.array(diff_vec)

# now multiply the diffusion matrix with this vector to get our ranked list of
# proteins
diff_result = diff_mat * diff_vec

print('Preparing output')
# doesn't come with labels but we know the order of the elements in diff_result
# and ppi_graph.nodes() is the same
full_df = pd.DataFrame(
    list(zip(ppi_graph.nodes(), diff_result)),
    columns = ['protein', 'rank']
)
full_df.to_csv('diffusion_result.csv', index = False)
