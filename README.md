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

## Model Building

## Web Access
