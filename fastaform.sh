
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END { printf("\n");}' < testfasta.rtf > tempfast.fasta
sed 's/[\]//g' tempfast.fasta > testfasta.fasta
