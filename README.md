# ZincBindPredict

[![Build Status](https://travis-ci.org/samirelanduk/ZincBindPredict.svg?branch=master)](https://travis-ci.org/samirelanduk/ZincBindPredict)
[![Coverage Status](https://coveralls.io/repos/github/samirelanduk/ZincBindPredict/badge.svg?branch=master)](https://coveralls.io/github/samirelanduk/ZincBindPredict
)

This is a machine learning project which predicts zinc binding residues from a
protein structure or protein sequence.

Its principal components are:

1. Scripts for generating training data from the ZincBindDB API. The ZincBind
database contains all known zinc binding structures, and these tools access the
database via its GraphQL API, and generate .csv files with the relevant
information.

2. Scripts for creating supervised learning models from the training data.

3. A small django application that provides access to the models over the web.

## Dataset Generation

The machine learning algorithms need training data so that they know how to
recognise a zinc binding site when they come across one. The actual object
that the classifiers will be given and asked to classify as either zinc binding
or not zinc binding is a set of residues.

Each family of zinc binding site will have its own classifier, and therefore its
own training data. So there will be the H3 classifer, which uses a H3.csv file
of training data to determine whether any three histidine residues it is given
could be a H3 binding site.

A combination of residues is represented by certain key metrics - distance
between them, etc. For each family there needs to be a structure training data
file and a sequence training data file.
## Model Building

The core of ZincBindPredict is the models - these are the things created from
the datasets created above, and the things accessed via the API below.

The models take as their input a set of residues and output true (zinc binding)
or false (not zinc binding). In practice a user would not submit a single
combination of residues, but rather a whole protein. The protein's potential
residue combinations would be identified, each one passed to the models in turn,
and then the protein's most likely combinations returned at the end.

The quality of the models will be measured by their recall and precision, and
the combination of these, the F1 Score. AUC/ROC scores will also be measured.

The datasets used to train the models will be those generated above, accessed as
.csv files using pandas.

## Web Access

While the user can download or clone this repository and generate their own
models, the 'official' ZincBind models will be hosted and provided to users,
accessible via an API. Thus this repository is also a django application.
