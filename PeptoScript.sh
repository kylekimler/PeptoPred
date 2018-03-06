echo -e "
Greetings and Welcome to PeptoPred. 
Use me to predict signal peptides, or maybe for some other binary classifications - as you like it"

sleep 2

read -p "whaddya wanna do? 
'1' to train a new model, '3' to predict unknowns: 
" theseareyouroptions

case $theseareyouroptions in 
	[1]* ) cd bin ; python3.6 PeptoModel.py ${trainbank} ${outputfile} ${windowsize}; cd .. ;;
	[2]* ) cd bin ; python3.6 PeptoCross.py ${setsize} ; cd .. ;; #look into using the same function for formatting for generating model as for testing. Consider writing cross training protocol in same py file
	[3]* ) cd bin ; python3.6 PeptoPred.py ${questionableproteins} ${outputfile} ; cd .. ;;
	[4]* ) echo -e "\e[1;40m" ; clear ; while :; do echo $LINES $COLUMNS $(( $RANDOM % $COLUMNS)) $(( $RANDOM % 72 )) ;sleep 0.05; done|awk '{ letters="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@#$%^&*()"; c=$4; letter=substr(letters,c,1);a[$3]=0;for (x in a) {o=a[x];a[x]=a[x]+1; printf "\033[%s;%sH\033[2;32m%s",o,x,letter; printf "\033[%s;%sH\033[1;37m%s\033[0;0H",a[x],x,letter;if (a[x] >= $1) { a[x]=0; } }}' ;;
	* ) echo "I guess you wanna do something else. Do you think this is Starbucks or something? There's no secret menu"
esac




