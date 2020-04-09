# network_diffusion.py
# take the COVID-19-human interactions, join them to a "normal" human PPI, and
# do network diffusion

import igraph
import scipy as sp
import numpy as np

# make a networkx graph object and add edges to it as we read the interactions
graph = igraph.Graph()

print('Adding human PPI edges')
# read in human PPI
i = 0
with open('string_human_ppi.txt', 'r') as in_file:
    # skip header
    x = next(in_file)
    for line in in_file:
        i += 1
        if i % 10000 == 0:
            print(f'On line {i}')
        cols = line.rstrip('\n').split(' ')
        # have to add vertices before corresponding edges with igraph
        graph.add_vertices(cols[:2])
        # returns the edge object it just added
        x = graph.add_edge(cols[0], cols[1], weight = cols[2])

# before reading in human-covid interactions, read in key for converting
# uniprot IDs to string IDs
id_dict = dict()
with open('string_ids_for_covid_interactors.csv', 'r') as in_file:
    for line in in_file:
        cols = line.rstrip('\n').split(',')
        # uniprot ID as key, string ID as value
        id_dict[cols[0]] = cols[1]

print('Adding human-COVID edges')
# now add in edges for human-covid interactions
# also make a list for later

# as written, this step takes about 6.5 hours. consider fixing that

interactors = list()
with open('supp_table_2.tsv', 'r') as in_file:
    # skip both header lines
    x = next(in_file)
    x = next(in_file)
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        covid_prot = cols[0]
        # convert IDs as we go; a few proteins had no STRING ID
        try:
            human_prot = id_dict[cols[7]]
            interactors.append(human_prot)
            score = cols[3]
            graph.add_vertices([human_prot, covid_prot])
            # returns copy of the edge just added
            x = graph.add_edge(covid_prot, human_prot, weight = score)
        except KeyError:
            continue

# compute diffusion matrix
print('Computing Laplacian Matrix')
lap_mat = graph.laplacian()
id_mat = sp.sparse.identity(lap_mat.shape[0])
print('Computing Diffusion Matrix')
diff_mat = np.linalg.inv(id_mat + (0.1 * lap_mat))

# need to create a binary vector with 1s for all COVID genes and 0s for all
# human genes in the same order as the rows in diff_mat
diff_vec = list()
# testing on small graphs indicates that the order of nodes in graph.nodes()
# always matches the order of columns in the Laplacian matrix
for prot in graph.nodes():
    if prot in interactors:
        # append lists with individual numbers so that the resulting list can
        # be coerced into a (1,n) numpy array
        diff_vec.append([1])
    else:
        diff_vec.append([0])

diff_vec = np.array(diff_vec)

diff_result = diff_vec * diff_mat
print(diff_result)
