"""Utility functions used in the generation of datasets."""

import random
import os
import subprocess
import re
import requests
from ftplib import FTP
from io import BytesIO
import gzip
import json
import kirjava
import atomium
from collections import Counter
import sys
sys.path.append("../zincbindpredict")

def parse_data_args(args, clustering=True):
    """Gets the desired settings from the args list."""

    similarity = [None]
    with open(os.path.join("data", "families.dat")) as f:
        families = f.read().splitlines()
    for arg in args:
        if arg.startswith("--limit="):
            families = [f for f in families if f in arg[8:].split(",")]
        if arg.startswith("--exclude="):
            families = [f for f in families if f not in arg[10:].split(",")]
        if arg.startswith("--clustering="):
            similarity = [None if s == "n" else float(s) for s in arg[13:].split(",")]
    return (families, similarity) if clustering else families


def fetch_data(url, query, variables):
    """Gets data from the ZincBindDB API and formats it to remove all the edges
    and nodes and stuff."""

    def strip_dict(d):
        for key, value in d.items():
            if isinstance(value, dict) and "edges" in value:
                d[key] = [strip_dict(e["node"]) for e in value["edges"]]
            elif isinstance(value, dict):
                d[key] = strip_dict(d[key])
        return d
    data = kirjava.execute(url, query, variables=variables)["data"]
    strip_dict(data)
    return list(data.values())[0]


def chars_to_family(chars):
    """Takes a list of characters and constructs a family from them. So, A1B2
    would be created from ['B', 'A', 'B'] for example."""

    counter = Counter(chars)
    return "".join(sorted([char + str(n) for char, n in counter.items()]))


def split_family(family):
    """Takes a family such as 'C3H1' and splits it into subfamilies such as 'C3'
    and 'H1'."""

    subfamilies, subfamily = [], ""
    for char in family:
        if char.isalpha() and subfamily:
            subfamilies.append([subfamily[0], int(subfamily[1:])])
            subfamily = ""
        subfamily += char
    subfamilies.append([subfamily[0], int(subfamily[1:])])
    return subfamilies


def save_csv(positives, negatives, name, similarity, path):
    """Takes a list of positive samples and a list of negative samples, and
    saves them to CSV."""

    if not positives and not negatives: return
    similarity = f"_{int(similarity * 100)}" if similarity else ""
    path = f"{path}{os.path.sep}{name}{similarity}.csv"
    lines = []
    flag = "w"
    lines.append(",".join((positives + negatives)[0].keys()) + ",positive")
    for index, samples in enumerate([positives, negatives]):
        for sample in samples:
            lines.append(",".join([str(v) for v in sample.values()] + [str(1 - index)]))
    with open(path, flag) as f:
        f.write("\n".join(lines) + "\n")


def cluster_sequences(sequences, similarity):
    """Takes a list of sequences and clusters them by some amount. CD-HIT must
    be installed."""

    lines = [f">{i}\n{s.upper()}" for i, s in enumerate(sequences)]
    try:
        with open("chains.fasta", "w") as f: f.write("\n".join(lines))
        word_size = 5 if similarity > 0.7 else 4 if similarity > 0.6 else\
            3 if similarity > 0.5 else 2
        subprocess.call(
            "cd-hit -i chains.fasta -d 0 -o temp -c {} -n {} -G 1 -g 1 -b 20 "
            "-s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 32000".format(similarity, word_size),
            shell=True, stdout=subprocess.PIPE
        )
        with open("temp.clstr") as f: clusters = f.read()
        clusters = clusters.split(">Cluster ")[1:]
        clusters = [[sequences[int(i)] for i in clust] for clust in reversed(
            sorted([
                re.compile(r">(.+?)\.\.\.").findall(c) for c in clusters
            ], key=len)
        )]
    finally:
        if os.path.exists("chains.fasta"): os.remove("chains.fasta")
        if os.path.exists("temp"): os.remove("temp")
        if os.path.exists("temp.clstr"): os.remove("temp.clstr")
    return [c[0] for c in clusters]


def sequence_contains_family(sequence, family):
    """Checks if a sequence contains a capitalied site of a given family."""

    residues = sorted([r for r in sequence if r.isupper()])
    return chars_to_family(residues) == family


def get_all_uniprot_sequences():
    """Gets all sequences un Uniprot, either from local disc (ideally), or from
    the Uniprot FTP server."""
    
    try:
        with open("uniprot_all.fasta") as f: text = f.read()
    except FileNotFoundError:
        ftp = FTP("ftp.uniprot.org")
        ftp.login()
        ftp.cwd("pub/databases/uniprot/current_release/knowledgebase/complete/")
        f = BytesIO()
        ftp.retrbinary("RETR uniprot_sprot.fasta.gz", f.write)
        f.seek(0)
        text = gzip.decompress(f.read()).decode()
        with open("uniprot_all.fasta", "w") as f: f.write(text)
    return ["".join(s.splitlines()[1:]).lower() for s in  text.split(">")[1:]]


def get_random_sequence_site(sequence, family):
    """Takes a sequence string and gets a random binding site for a given
    family. If there is no matching site, returns None."""

    indices = []
    for subfamily in split_family(family):
        code, count = subfamily[0].lower(), subfamily[1]
        subfamily_indices = [index for index, char in
            enumerate(sequence.lower()) if char == code]
        if len(subfamily_indices) < count: return None
        indices += random.sample(subfamily_indices, count)
    return "".join([char.upper() if index in indices else char.lower()
        for index, char in enumerate(sequence)])


def get_atomium_site(site):
    """Takes the API JSON for a site and gets a list of atomium residues
    representing that site."""

    pdb = atomium.fetch(site["pdb"]["id"])
    model = pdb.generate_assembly(site["pdb"]["assembly"])
    model.optimise_distances()
    site_residues = [r for r in site["residues"] if r["chainSignature"]]
    return get_residues_from_model(model, site_residues)


def get_residues_from_model(model, residues):
    """Takes a model, and some residues in JSON form, and then finds the
    matching residues in the model. Sometimes this function fails."""

    at_residues = []
    for res in residues:
        possibles = model.residues(id=res["atomiumId"])
        matching = [r for r in possibles if r.atom(name="CA").location == (
            res["atoms"][0]["x"], res["atoms"][0]["y"], res["atoms"][0]["z"]
        )]
        at_residues.append(matching[0])
    return at_residues


def get_all_pdb_codes():
    """Get all PDB codes from the RCSB."""

    zinc_codes = [p["id"] for p in fetch_data(
        "https://api.zincbind.net", "{ pdbs { edges { node { id } } } }", {}
    )]
    query = {
        "query": {"type": "terminal", "service": "text"},
        "request_options": {"pager": {"start": 0, "rows": 1000000000}},
        "return_type": "entry"
    }
    url = "https://search.rcsb.org/rcsbsearch/v1/query?json="
    return [pdb["identifier"] for pdb in requests.get(
        url + json.dumps(query)
    ).json()["result_set"] if pdb["identifier"] not in zinc_codes]


def get_random_atomium_site(code, family):
    """Gets a random binding site in the form of atomium residues matching some
    family. If the code given has no matching site, None is returned."""

    pdb = atomium.fetch(code)
    model = atomium.fetch(code).model
    model.optimise_distances()
    residues = []
    for subfamily in split_family(family):
        code, count = subfamily[0].upper(), subfamily[1]
        matching_residues = list(model.residues(code=code))
        if len(matching_residues) < count: return None
        residues += random.sample(matching_residues, count)
    return residues