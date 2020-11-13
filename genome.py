import joblib
import os
import time
from tqdm import tqdm
import warnings
from server.utilities import *
from data.common import sequence_site_to_vector
warnings.warn = lambda *args, **kwargs: None

from random import randint
import requests
import json
from ftplib import FTP
from io import BytesIO
import gzip

URL = ("https://www.ncbi.nlm.nih.gov/Structure/ngram?&callback=jQuery22407272655"
"846776386_1604411050391&start={}&limit={}&q=%5Bdisplay()%2Chist(kingdom%2Cgrou"
"p%2Csubgroup%2Clevel%2Cpartial%2Canomalous%2Chost%2Crefseq_category)%5D.from(G"
"enomeAssemblies).usingschema(%2Fschema%2FGenomeAssemblies).matching(tab%3D%3D%"
"5B%22Prokaryotes%22%5D).sort(organism%2Casc)&_search=false&rows=50&page=3&sidx"
"=organism&sord=asc&_=1604411050395")

models = {
    f.split("-")[0]: joblib.load("predict/models/sequence/" + f)
    for f in os.listdir("predict/models/sequence") if f.endswith("joblib")
}

def check_proteins(proteins):
    zbps = []
    for protein in tqdm(proteins):
        if len(protein) < 1000:
            for family, model in models.items():
                try:
                    possibles = list(sequence_to_family_inputs(protein, family))
                    dicts = [sequence_site_to_vector(possible) for possible in possibles]
                    vectors = [list(d.values()) for d in dicts]
                    if not vectors: continue
                    max_length = max(len(v) for v in vectors)
                    vectors = [v for v in vectors if len(v) == max_length]
                    probabilities = [p[1] for p in model.predict_proba(vectors)]
                    if any(p >= 0.99 for p in probabilities):
                        zbps.append(protein)
                        break
                except Exception: pass
    return len(zbps) / len(proteins)


def print_and_save(s):
    with open("output.log") as f:
        lines = [l + "\n" for l in f.read().splitlines()]
        lines.append(s + "\n")
        with open("output.log", "w") as f:
            f.write("".join(lines))
            print(s)
    



while True:
    integer = randint(1, 100000)
    data = json.loads(requests.get(URL.format(integer, 1)).content[41:-2])
    genome = data["ngout"]["data"]["content"][0]

    if "ftp_path_refseq" in genome:
        ftp = FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        path = genome["ftp_path_refseq"][27:]
        ftp.cwd(path)
        files = ftp.mlsd(".")
        for file_name, _ in files:
            if "translated_cds.faa" in file_name:
                f = BytesIO()
                ftp.retrbinary("RETR " + file_name, f.write)
                f.seek(0)
                bytestring = f.read()
                contents = gzip.decompress(bytestring).decode()
                genes = contents.split(">")[1:]
                proteins = ["".join(g.splitlines()[1:]) for g in genes]
                if len(proteins) < 4000:
                    print_and_save(genome["organism"] + ", " + str(genome["id"]))
                    ratio = check_proteins(proteins[:int(len(proteins) / 10)])
                    print_and_save(f"{round(ratio * 100, 1)}%")
                    if ratio < 0.14:
                        ratio = check_proteins(proteins)
                        print_and_save(f"{round(ratio * 100, 1)}%")
                    print_and_save("")
                    break
        
                
            
            
        '''with open('README', 'wb') as fp:
            ftp.retrbinary(path, fp.write)'''


        


'''genomes = [f for f in os.listdir() if f.endswith(".faa")]

for genome in genomes:
    print()
    print(genome)

    with open(genome) as f:
        genome = f.read()

    genes = genome.split(">")[1:]
    proteins = ["".join(g.splitlines()[1:]) for g in genes]
    print("There are", len(proteins), "proteins")

    models = {
        f.split("-")[0]: joblib.load("predict/models/sequence/" + f)
        for f in os.listdir("predict/models/sequence") if f.endswith("joblib")
    }

    zbps = []
    for protein in tqdm(proteins[:5]):
        
        for family, model in models.items():
            possibles = list(sequence_to_family_inputs(protein, family))
            dicts = [sequence_site_to_vector(possible) for possible in possibles]
            vectors = [list(d.values()) for d in dicts]
            if not vectors: continue
            max_length = max(len(v) for v in vectors)
            vectors = [v for v in vectors if len(v) == max_length]
            probabilities = [p[1] for p in model.predict_proba(vectors)]
            if any(p == 1 for p in probabilities):
                zbps.append(protein)
                break
    print(len(zbps), len(proteins), round(len(zbps) / len(proteins) * 100))'''


