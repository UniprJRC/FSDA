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

for fold in clustering combinatorial examples FSDAdemos graphics multivariate regression robust utilities
do
	cd $fold
	>patrizia
	mv Contents.m patrizia
	echo "% $fold" | tr '[:lower:]' '[:upper:]' >Contents.m
	echo "%" >>Contents.m
	echo "% Files :" >>Contents.m

	ls -1 *.m |grep -v "Contents.m" | while read file
	do
		riga=`grep -m 1 "^%" $file`
		comm=`echo "$riga" |cut -c2-`
		base=`basename $file .m`
		printf "%% %-28s - %s\n" "$base" "$comm" >>Contents.m
	done
	unix2dos Contents.m
	
	echo "Contents.m in $fold GENERATO !!!!"
	echo "DIFFERENZA  di Contents.m < nuovo_file > vecchio_file"
	sdiff Contents.m patrizia
	echo " <Enter> per continuare "
	read avanti
	cd ..
done

cd datasets
>patrizia
mv Contents.m patrizia

echo "<a href=""matlab: docsearchFS('datasets')"">Link to the help function</a>" >Contents.m
echo "<a href=""matlab: docsearchFS('datasets_mv')"">Link to the help function</a>" >>Contents.m
echo "" >>Contents.m

find clustering multivariate multivariate_regression regression -name "*.txt" >>Contents.m
sed -e ' {
        1,2 s/^/%/
        3,$ s/^/% /
        }' Contents.m
unix2dos Contents.m

echo "Contents.m in datasets GENERATO !!!!"
echo "DIFFERENZA  di Contents.m < nuovo_file > vecchio_file"
sdiff Contents.m patrizia
echo " <Enter> per continuare "
read avanti
echo " I file Contents.m sono stati ricreati. I vecchi Contents.m sono salvati come 'patrizia'"

exit 0

