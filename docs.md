# Server

The ZincBindPredict interface is a GraphQL API. You perform searches by sending either a searchSequence mutation, or a searchStructure mutation. In return you get a job ID, which can be used to query the job as it progresses.

## searchSequence

The searchSequence mutation takes a raw (that is - not FASTA, no line breaks etc.) protein sequence to be searched as the sequence argument. It also takes optional arguments of families to filter by. It returns the job ID that can be used to query the sequence job.

For each family (either all of them or the ones provided), all possible combinations of residues matching that family will be identified, these combinations will be turned into a feature vector relevant to the family's sequence model, and each of these potential sites will be passed to that model one by one.

As the job progresses, a list of predicted sites and rejected sites will be built up, with each site having a probability, family, and a representation of the residues in question.

The sequence job can be queried at any time - it will give the status of the sequence job (which family it is looking through currently), the time it was started, a string representation of the protein it was given, and a list of predicted and rejected sequences in the form outlined in the previous paragraph.

## searchStructure

The searchStructure mutation takes a protein structure to be searched, and optionally a list of families to filter by, like the searchSequence mutation. Unlike that one however, the structure is supplied as a file upload, not a string. In addition, it contains three boolean flags - useFamiliesModels, useLocationModels, and useHalf. By default structure searching looks for whole binding sites using family based models, whole binding sites using location based models, and half binding sites using location based models. These flags turn off these searches respectively.

Assuming all of these are turned on, the structure searching job will proceed as follows.

...continue...
