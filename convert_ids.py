# convert_ids.py

# STRING used Ensembl protein IDs and the preprint used Uniprot IDs

# Uniprot IDs to multiple Ensembl protein IDs, so we're converting all Ensembl
# protein IDs to Uniprot IDs (turns out all of the Ensembl protein IDs we have
# from STRING map onto unique uniprot IDs, which is great)

# There are tons of IDs so we're doing this programmatically

import urllib.parse
import urllib.request

# first read in the STRING file and turn it into a list of all the unique IDs
# we want to convert
# use a set to avoid duplicates
ensembl_ids = set()
with open('string_human_ppi.txt') as in_file:
    for line in in_file:
        # strip newlines and split on space
        cols = line.rstrip('\n').split(' ')
        # ignore first line
        if cols[0] == 'protein1':
            continue
        else:
            # third column is interaction strength
            for ensembl_id in cols[:2]:
                # all ensembl IDs are prefixed by 9606. for some reason
                ensembl_ids.add(ensembl_id.lstrip('9606.'))

# make the set into a delimited string
id_string = ' '.join(ensembl_ids)

# now interact with the uniprot API
url = 'https://www.uniprot.org/uploadlists/'

params = {
    # ensembl protein ID
    'from': 'ENSEMBL_PRO_ID',
    # normal uniprot ID
    'to': 'ID',
    'format': 'tab',
    'query': id_string
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8'))
