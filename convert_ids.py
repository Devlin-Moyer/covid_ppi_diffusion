# convert_ids.py

# We have STRING IDs for the human PPI and Uniprot IDs so we're converting
# the STRING IDs to Uniprot IDs. There's a suspiciously large number of unique
# STRING IDs (i.e. more than the number of unique proteins in the human genome)
# so most don't map to uniprot IDs, but that's fine

import urllib.parse
import urllib.request

# first read in the STRING file and turn it into a list of all the unique IDs
# we want to convert
# use a set to avoid duplicates
string_ids = set()
with open('string_human_ppi.txt') as in_file:
    for line in in_file:
        # strip newlines and split on space
        cols = line.rstrip('\n').split(' ')
        # ignore first line
        if cols[0] == 'protein1':
            continue
        else:
            # third column is interaction strength
            string_ids.add(cols[0])
            string_ids.add(cols[1])

# make the set into a space-delimited string
id_string = ' '.join(string_ids)

# now interact with the uniprot API
url = 'https://www.uniprot.org/uploadlists/'

params = {
    # ensembl protein ID
    'from': 'STRING_ID',
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
