#!/bin/bash
#set -x

#Installation script FSDA Toolbox for LINUX system


echo "+++++++++++++++== FSDA Toolbox Installation on Linux system ==+++++++++++++++++"
echo ""
echo "License and copyright information:"
echo ""
echo "The FSDA toolbox is protected by a European Union Public Licence (EUPL)."
echo "The EUPL is the first European Free/Open Source Software (F/OSS) Licence."
echo "For more information on the EUPL, please visit http://www.osor.eu/eupl"
echo ""
echo "Requisite for FSDA toolbox is MATLAB, which is copyright by © 1994-2017 The MathWorks, Inc."
echo ""
echo "AUTHORS and DEVELOPERS:"
echo ""
echo "	University of Parma: Marco Riani"
echo "						 Aldo Corbellini"
echo "						 Andrea Cerioli"
echo "						 Gianluca Morelli"
echo "						 Fabrizio Laurini"
echo ""						 
echo "	European Commission, Joint Research Centre: Domenico Perrotta"
echo "												Francesca Torti"
echo "												Emmanuele Sordini"
echo "												Patrizia Calcaterra"
echo "												Andrea Cerasa"
echo "												Massimiliano Gusmini"
echo ""
echo "  University of Milan: Silvia Salini (from 2016)"
echo ""                    	
echo "The sources are open and the use is regulated by EUPL."
echo ""


CURDIR=`pwd`


echo " Where do you want to install FSDA Toolbox ?"
echo " Please insert the complete path (default: /usr/bin/FSDA) --> "
read where

if [ "$where" == "" ]
then
	where="/usr/bin/FSDA"
else
	where="$where/FSDA"
fi


if [ -e $where ]
then 
	risp="yes"
	echo $where " already exists. Do you want to continue ? (yes/no) (default : yes) --> "
	read risp
	if [ $risp == "no" -o $risp == "No" -o $risp == "NO" -o $risp == "N" -o $risp == "n" ]
	then
		echo "FSDA Toolbox installation aborted ++++++++"
		exit 1
	fi
else
	mkdir -p $where
fi


cd $CURDIR/FSDA

echo "Copying FSDA Toolbox files in : " $where "............"

cp -R * $where/.


matlab -logfile $where/log.txt -nodesktop -nosplash -r "a=ver('MATLAB'); disp([a.Version]); quit;"


VER=`sed -n '$p' $where/log.txt`
REL=`echo $VER | cut -c1`

MATEXE=`which matlab`
MATPAT=${MATEXE%%/bin/matlab}

# matlab -wait -automation -nodesktop -r " builddocsearchdb ('${where}/FSDA/helpfiles/FSDA') ; quit"

if [ $REL -ge 8 ] 
then
    mv $where/helpfiles/FSDA $where/helpfiles/FSDAtomove
	mkdir $where/helpfiles/FSDA
	mv $where/helpfiles/FSDAtomove/helptoc.xml $where/helpfiles/FSDA/.
	mv $where/helpfiles/FSDAtomove/fsda_product_page.html $where/helpfiles/FSDA/.
	mv $where/helpfiles/FSDAtomove ${MATPAT}/help/.
	mv ${MATPAT}/help/FSDAtomove ${MATPAT}/help/FSDA
#	matlab -nodesktop -r " addpath ${MATPAT}/help/FSDA ; builddocsearchdb ('${MATPAT}/help/FSDA') "
fi


echo "Setting MATLAB environment ..."	

matlab -nodesktop -nojvm -r "addpath '$where/examples' ; addpath '$where/utilities' ; addpath '$where/combinatorial' ; addpath '$where/FSDAdemos' ; addpath '$where/graphics' ; addpath '$where/utilities_stat' ; addpath '$where/utilities_help' ; addpath '$where/datasets/multivariate' ; addpath '$where/datasets/regression' ; addpath '$where/datasets/multivariate_regression' ; addpath '$where/datasets/clustering' ; addpath '$where/clustering' ; addpath '$where/regression' ; addpath '$where/multivariate' ; addpath '$where' ; savepath ; exit "


if [ $REL -ge 8 ] 
then
	matlab -nodesktop -r " cd $where; matlab.apputil.install('brushRES'); quit;"
	matlab -nodesktop -r " cd $where; matlab.apputil.install('brushFAN'); quit;"
	matlab -nodesktop -r " cd $where; matlab.apputil.install('brushROB'); quit; "
fi 

echo ""
echo ""
echo "The following folders are added to MATLAB path :"
echo ""
echo "{FSDA toolbox}"
echo "{FSDA toolbox}/clustering"
echo "{FSDA toolbox}/regression"
echo "{FSDA toolbox}/multivariate"
echo "{FSDA toolbox}/utilities_stat"
echo "{FSDA toolbox}/utilities_help"
echo "{FSDA toolbox}/FSDAdemos"
echo "{FSDA toolbox}/graphics"
echo "{FSDA toolbox}/utilities"
echo "{FSDA toolbox}/examples"
echo "{FSDA toolbox}/combinatorial"
echo "{FSDA toolbox}/datasets/clustering"
echo "{FSDA toolbox}/datasets/multivariate"
echo "{FSDA toolbox}/datasets/regression"
echo "{FSDA toolbox}/datasets/multivariate_regression"
echo ""
echo "Please do not remove them otherwise the FSDA toolbox will not work."
echo ""
echo "For more information, please read $where/InstallationNotes.pdf"
echo " +++++++++++++++++++== FSDA Toolbox Installation end ==++++++++++++++++++"
