#!/bin/env python3

import urllib.parse
import urllib.request
#import Bio.SeqIO
import subprocess
import argparse
import os.path
import time

parser = argparse.ArgumentParser()
parser.add_argument('file', help='A required file with UPIs')
args = parser.parse_args()

# get the UPIs from file
inp = args.file #'UPIs/upi_00'
with open(inp, 'r') as input_handle:
    upis = [line.strip() for line in input_handle.readlines()]

upi_str = " ".join(upis)

outp = os.path.join('filtered_results','uniprot','seq_' + os.path.basename(inp)[4:])

url = 'https://www.uniprot.org/uploadlists/'

# map the UPIs to Uniprot entry names
params = {
'from': 'UPARC',
'to': 'ACC',
'format': 'tab',
'query': upi_str #'P40925 P40926 O43175 Q9UM73 P97793'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
#with urllib.request.urlopen(req) as f:
#    response = f.read()

try:
    with urllib.request.urlopen(req) as f:
        response = f.read()
except:
   time.sleep(60)
   with urllib.request.urlopen(req) as f:
       response = f.read() 


mapping = response.decode('utf-8')
ids = [line.split('\t') for line in mapping.split('\n')]
ids = [id_entry[1]  for id_entry in ids if len(id_entry) == 2][1:]


if len(ids) == 0:
    subprocess.run(['touch',outp])
elif len(ids) == 1:
    url = 'uniprot.org/uniprot/?format=tab&columns=id,organism,sequence&query=accession%3A'+ids[0]
    subprocess.run(['wget', '-O', outp, url])
else:
    url = 'uniprot.org/uniprot/?format=tab&columns=id,organism,sequence&query=accession%3A'+ids[0]
    for id_entry in ids[1:]:
        url += '+OR+accession%3A'+id_entry

    subprocess.run(['wget', '-O', outp, url])
  

# retrive corresponing organism and sequence
# 'uniprot.org/uniprot/?format=tab&columns=id,organism,sequence&query=accession%3A{ID1}+OR+accession%3A{ID2}+OR+accession%3A{ID3}'
#url = 'uniprot.org/uniprot/?format=tab&columns=id,organism,sequence&query=accession%3A'+ids[0]
#for id_entry in ids[1:]:
#    url += '+OR+accession%3A'+id_entry

#subprocess.run(['wget', '-O', outp, url])


