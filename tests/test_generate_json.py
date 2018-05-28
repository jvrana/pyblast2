import json
import requests
import os
from glob import glob
from tqdm import tqdm

query = """
query ParseSequences($file: String!, $filename: String!) {
  tojson(file: $file, filename: $filename) {
    messages
    parsedSequence {
      name
      sequence
      circular
      description
      features {
        name
        type
        id
        start
        end
        strand
      }
    }
    success
  }
}
"""

def test_generate_subjects(here):

    def parse_to_json(file, filename):
        url = 'http://0.0.0.0:3000/graphql'
        variables = {
            "file": file,
            "filename": filename
        }
        payload = {"query": query, "variables": variables}
        headers = {'content-type': 'application/json'}
        return requests.post(url, data=json.dumps(payload), headers=headers)


    # generate templates
    directory = os.path.join(here, "data/test_data/genbank/templates")

    sequences = []
    for path in tqdm(glob(os.path.join(directory, "*.gb"))):
        with open(path, 'r') as f:
            res = parse_to_json(f.read(), os.path.basename(path))
            sequences.append(res.json())

    sequences_json = [seq["data"]["tojson"]["parsedSequence"] for seq in sequences]
    with open(os.path.join(here, "data/test_data/templates.json"), 'w') as f:
        json.dump(sequences_json, f)

    # generate query
    with open(os.path.join(here, "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb"), "r") as f:
        res = parse_to_json(f.read(), f.name)
        with open(os.path.join(here, "data/test_data/query.json"), 'w') as out:
            json.dump(res.json()["data"]["tojson"]["parsedSequence"], out)
