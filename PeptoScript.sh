echo -e "
Greetings and Welcome to PeptoPred. 
PeptoPred comes with a training dataset bestowed upon us by the left hand of God
It also comes with a trained model made by Random Forest Classification and an amino acid window size of 15
If you want to make a model, run $python3.6 PeptoModel.py {x} where x is the desired window size
PeptoModel.py is currently set up to create Decision Tree Classifier models, but you can edit it to create the model of your desired

RandomForest modeling of the whole included dataset takes about 14 minutes on a macbook air

All created models will appear in the PeptoPred/Models folder.
----------------------------------------------------------------------

If you'd like to test some peptides on the given model, 
input the name of a fastafile and this program will output a prediction into the file 'results/[yourfilename]RFC15Pred.txt'

You can easily change the PeptoPred.py file to i/o with your new model - lines 12 and 46

Otherwise, "

read fastain

python3.6 PeptoPred.py ${fastain}

