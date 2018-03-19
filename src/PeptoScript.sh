echo -e "
Greetings and Welcome to PeptoPred. 
Right now we have a RandomForestClassifier model of the whole dataset, with a window size of 15. 
RandomForest modeling of the whole dataset takes about 14 minutes on a macbook air
-----------
In the future, expect to see new feature extraction methods for better training, the option to train on your own dataset, and a cross training option.

For now,

If you'd like to test some peptides on the model, input the name of a fastafile and this program will output a prediction into the file '[yourfilename]Pred.txt'"

read -p "whaddya wanna do? 
'1' to train a new model on the included dataset, '2' to cross train a model, '3' to predict unknowns: 
" theseareyouroptions

case $theseareyouroptions in 
	[1]* ) python3.6 PeptoModel.py '../data/globular_signal_peptide_2state.3line.txt' 15 ;;
	[2]* ) python3.6 PeptoCross.py ${setsize} ;; 
	[3]* ) python3.6 PeptoPred.py ${questionableproteins} ${outputfile} ;;
	[4]* ) clear ; while :; do echo $LINES $COLUMNS $(( $RANDOM % $COLUMNS )) $(( $RANDOM % 72 )) ;sleep 0.05; done|awk '{ letters="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@#$%^&*()"; c=$4; letter=substr(letters,c,1);a[$3]=0;for (x in a) {o=a[x];a[x]=a[x]+1; printf "\033[%s;%sH\033[2;32m%s",o,x,letter; printf "\033[%s;%sH\033[1;37m%s\033[0;0H",a[x],x,letter;if (a[x] >= $1) { a[x]=0; } }}' ;;
	* ) echo "I guess you wanna do something else. Do you think this is Starbucks or something? There's no secret menu"
esac

