#!/bin/sh

#set -x
echo "Inserisci il nome dell'archivio FSDA da analizzare e aggiornare : "
read what

if [ ! -d "$what" ]
then
        echo "ERRORE :" $what " non esiste !!!!"
        exit 1
fi

cd $what

for fold in clustering combinatorial examples graphics multivariate regression robust utilities
do
	cd $fold
	echo "% $fold" | tr '[:lower:]' '[:upper:]' >Contents.m
	echo "%" >>Contents.m
	echo "% Files :" >>Contents.m

	ls -1 *.m |grep -v "Contents.m" | while read file
	do
		riga=`grep -m 1 "^%" $file`
		comm=`echo $riga |cut -c2-`
		base=`basename $file .m`
		printf "%% %-28s - %s\n" "$base" "$comm" >>Contents.m
	done
	unix2dos Contents.m
#echo "DIFFERENZA PER $fold"
#sdiff Contents.m patri.m
echo Contents.m in $fold generato
echo prosegui
read avanti
	cd ..
done

cd datasets
echo prosegui
read avanti

echo "<a href=""matlab: docsearch('datasets_reg')"">Link to the help function</a>" >Contents.m
echo "<a href=""matlab: docsearch('datasets_mult')"">Link to the help function</a>" >>Contents.m
echo "" >>Contents.m

find multivariate multivariate_regression regression -name "*.txt" >>Contents.m
sed -e ' {
        1,2 s/^/%/
        3,$ s/^/% /
        }' Contents.m
unix2dos Contents.m
echo Contents.m in datasets generato

exit 0

