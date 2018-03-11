echo -e "
Greetings and Welcome to PeptoPred. 
Right now we have a RandomForestClassifier model of the whole dataset, with a window size of 15. 
RandomForest modeling of the whole dataset takes 764 seconds on 8gb RAM. 
-----------
In the future, expect to see new feature extraction methods for better training, the option to train on your own dataset, and a cross training option.

For now,

If you'd like to test some peptides on the model, input the name of a fastafile and this program will output a prediction into the file '[yourfilename]Pred.txt'"

read fastain

python3.6 PeptoPred.py ${fastain}

