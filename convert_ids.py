# convert_ids.py

# get the string IDs for all of the human proteins found in the preprint

import requests

# get list of uniprot IDs for the human proteins that interacted with COVID
# proteins

prot_list = list()
with open('supp_table_2.tsv', 'r') as in_file:
    # skip both header rows
    x = next(in_file)
    x = next(in_file)
    for line in in_file:
        uniprot_id = line.rstrip('\n').split('\t')[7]
        prot_list.append(uniprot_id)

# interact with STRING API

# construct URL

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

request_url = "/".join([string_api_url, output_format, method])

params = {
    "identifiers" : '\r'.join(prot_list),
    "species" : 9606, # human NCBI identifier 
    "limit" : 1, # only one (best) identifier per input protein
    "echo_query" : 1, # see input identifiers in output
    "caller_identity" : "dcmoyer@bu.edu" # they want a personal identifier
}

results = requests.post(request_url, data=params)

# write output
with open('string_ids_for_covid_interactors.csv', 'w') as out:
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        input_identifier, string_identifier = l[0], l[2]
        out.write(input_identifier + ',' + string_identifier + '\n')
