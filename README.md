# Repeat-analysis-pipeline
A complete software pipeline for analyzing repeat content in genomes.

Repeat-analysis-pipeline is a complete repeat pipeline which detects repeat elements in a given genome and also classifies the repeat elements into different classes.  RepeatModeler is used for de novo detection of repeat elements. This uses two de novo repeat searching algorithms RECON and RepeatScout. It then uses CENSOR for annotating these repeat elements using the Repbase library as reference. Interproscan is also used for searching specific protein domains which may be specific to a certain type of transposable element.

Once the repeat library is created, it is run using repeatmasker throughout the entire genome to determine the locations of each of these repeat elements. 

A second round of annotations of the repeat elements is carried out using a machine learning type approach. It searches for pentanucleotide frequencies (4*4*4*4*4=1024 unique combinations) in each of the DNA sequences of the repeat elements. The frequencies of each of these repeat elements is then used as a signature for classifying unknown repeat elements. This algorithm follows a decision tree approach where in each step the algorithm a repeat element gets classified according to the hierarchy of repeat elements (such as gypsy is a subset of LTR retrotransposons which in turn is a subset of retrotransposons). The algorithm uses grid search for cross-validating the training dataset and also reports the confusion matrix and accuracy of the classification using a validation step.
