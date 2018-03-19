# PeptoPred

To run PeptoPred on a fasta of your choice, change directory to src and run "PeptoPred.py [fastafilename]". You can use 2 or 3 line fastas on it interchangeably. There's a stats comparison piece to the program you can activate if you want to see accuracy for yourself, otherwise there's some stats in the report. 

There are a couple test fastas in the data folder - "data/testpeps.txt" and "data/testpeps2.txt" are 3-liners that come from the given dataset, while testfasta.fasta is the beginning of swissprot (2 line) reformatted with "fastaform.sh" to fit the feature extractor "PeptoPred/src/Formatting.py". 

To run on testfasta.fasta - go to the src folder and paste "python3.6 PeptoPred.py ..data/testfasta/fasta" into the terminal

PeptoPred.py will print results to "results/[fastafilename][windowsize]Pred.txt"

You can run PeptoModel.py [fastafilename] [windowsize] to train a new model. You also have to be in the src directory for it to run properly. It will work on the dataset "../data/globular_signal_peptide_2state.3line.txt". It's currently set to create a decision tree classifer - you can change that to a different algorithm by changing line 36. 

You can also run PeptoCross.py [windowsize] to cross validate a decision tree classifier model - you also have to be in the src directory for that to work and you need to have previously created a DTC model with the chosen windowsize. It will print results to the terminal.


