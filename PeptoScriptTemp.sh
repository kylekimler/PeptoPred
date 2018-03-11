case $theseareyouroptions in 
	[1]* ) python3.6 PeptoModel.py ${trainbank} ${outputfile} ${windowsize} ;;
	[2]* ) python3.6 PeptoCross.py ${setsize} ;; #look into using the same function for formatting for generating model as for testing. Consider writing cross training protocol in same py file
	[3]* ) python3.6 PeptoPred.py ${questionableproteins} ${outputfile} ;;
	[4]* ) echo -e "\e[1;40m" ; clear ; while :; do echo $LINES $COLUMNS $(( $RANDOM % $COLUMNS)) $(( $RANDOM % 72 )) ;sleep 0.05; done|awk '{ letters="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@#$%^&*()"; c=$4; letter=substr(letters,c,1);a[$3]=0;for (x in a) {o=a[x];a[x]=a[x]+1; printf "\033[%s;%sH\033[2;32m%s",o,x,letter; printf "\033[%s;%sH\033[1;37m%s\033[0;0H",a[x],x,letter;if (a[x] >= $1) { a[x]=0; } }}' ;;
	* ) echo
esac




