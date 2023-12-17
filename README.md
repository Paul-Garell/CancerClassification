# CancerClassification

*Created by:*
Paul Garell (prg74)
Pedro Vera Perez (ppv7)
Ryan Mistretta (rhm255)


##Background

This repository contains a breast cancer type classification algorithm and model. Our function takes as input a sequence of DNA pairs. Internally this sequence is aligned over a large reference genome (~400,000 pairs). Once we find the most likely alignment, a series of single nucleotide polymorphisms are generated. With the series of SNPs, our trained model performs a series of inferences to determine the classification of cancer. 

Since we match the input DNA sequence with a fairly large reference genome, running the tests takes a long time (about 60 minutes on M2 mac). 

We provide both a trained model as well as our training script. The training script is highly modular and can be adapted to train on different types of cancer with ease. Additionally, we provide our sequence alignment code and methods to perform inferences on our model. The file inference.py contains functions for inferring on the model. Additionally the file E2Etest.py contains functions for performing a full pipeline inference. Currently it is set up to read from a text file in the format ["sequence", "Cancer type"]. We provide two test suites. Our full 41 test suite, and a reduced 10 sequence suite. The 10 sequence suite achieves a score of 7/10, and our full test suite achieves 36/41. 

##Runbook

To run our test suite: run 


