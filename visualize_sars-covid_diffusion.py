# visualize_sars-covid_diffusion.py

import networkx as nx
import pygraphviz as gv
import re

graph = gv.AGraph(splines = 'true')

print('Creating human-COVID-SARS PPI')

# start by getting all human-SARS interactions
# also make a list of all of these so we can separate their diffusion scores
# from the other proteins' scores (actually a set so elements are unique)
interactors = set()
with open('sars_human_interactions.tsv', 'r') as in_file:
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        # uniprot IDs are prefixed by 'uniprotkb:' 
        prot1 = cols[0].lstrip('uniprotkb:')
        prot2 = cols[1].lstrip('uniprotkb:')
        interactors.add(prot1)
        interactors.add(prot2)
        # the human and SARS proteins aren't segregated into different columns
        # so have to figure out which is which to color nodes properly
        if prot1.endswith('HUMAN'):
            graph.add_node(prot1, color = 'blue')
            # prot2 is the SARS protein
            graph.add_node(prot2, color = 'orange')
        else:
            # prot1 is the SARS protein
            graph.add_node(prot1, color = 'orange')
            graph.add_node(prot2, color = 'blue')
        graph.add_edge(prot1, prot2)

# add human-COVID interactions
with open('supp_table_2.tsv', 'r') as in_file:
    # skip both header lines
    x = next(in_file)
    x = next(in_file)
    for line in in_file:
        cols = line.rstrip('\n').split('\t')
        # nobody likes whitespace
        covid_prot = re.sub(' ', '_', cols[0])
        human_prot = cols[7]
        interactors.add(human_prot)
        interactors.add(covid_prot)
        graph.add_node(covid_prot, color = 'green')
        graph.add_node(human_prot, color = 'blue')
        graph.add_edge(covid_prot, human_prot)

print('Working with protein diffusion scores')
# read in diffusion results so we can pick out which of the many proteins in
# the full human PPI have particularly high scores
diffusion_scores = list()
with open('sars_diff.csv', 'r') as in_file:
    # skip header
    x = next(in_file)
    for line in in_file:
        diffusion_scores.append(line.rstrip('\n').split(','))

# scores for all proteins that aren't COVID proteins or their direct
# interaction partners
non_int_scores = [
    score for score in diffusion_scores if score[0] not in interactors
]
# sort by diffusion scores (each sublist is protein_name, diffusion_score)
non_int_scores_sorted = [
    sorted(sublist, key = lambda x: x[1]) for sublist in non_int_scores
]
# now just pick a smallish number of these proteins (second element of sublist)
top_other_proteins = [sublist[1] for sublist in non_int_scores_sorted[0:25]]

print('Adding edges for high-scoring human proteins')
# now get all of those proteins' interactions
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
        # skip edges with low scores
        if int(cols[2]) < 700:
            continue
        # convert IDs, but since most IDs won't map, use a try/except loop
        try:
            uniprot_ids = [id_dict[cols[0]], id_dict[cols[1]]]
            # only add the edge if one of the nodes is in top_other_proteins
            # and the other is a direct interaction partner of a COVID protein
            if ((uniprot_ids[0] in top_other_proteins and 
                uniprot_ids[1] in interactors) or
                (uniprot_ids[1] in top_other_proteins and
                uniprot_ids[0] in interactors)):
                # highlight the protein that is not an interactor
                if uniprot_ids[0] in top_other_proteins:
                    graph.add_node(uniprot_ids[0], color = 'red')
                elif uniprot_ids[1] in top_other_proteins:
                    graph.add_node(uniprot_ids[1], color = 'red')
                graph.add_edge(
                    uniprot_ids[0],
                    uniprot_ids[1]
                )
        except KeyError:
            continue
 
# Now drop all nodes with a degree of 1
for (node, degree) in graph.degree_iter():
    if degree == 1:
        graph.delete_node(node)

print('Creating image of network')
graph.draw('sars-covid_diff.png', prog = 'fdp')
