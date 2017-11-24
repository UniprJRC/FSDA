rem @echo off
set where=%1

set mathpath=%2

set drv=%where:~1,2%

%drv%

cd "%where%"

rem echo "SONO IN mgmhlpR8 . I parametri sono  " %where% " e " %mathpath%

cd FSDA
mkdir _tmp_helpfiles
cd helpfiles
rename FSDA FSDAtomove
mkdir FSDA
cd FSDAtomove
move helptoc.xml ..\FSDA
move fsda_product_page.html ..\FSDA

cd ..
move FSDAtomove ..\_tmp_helpfiles
cd ..\_tmp_helpfiles
rename FSDAtomove FSDA

rem echo "dovresti avere FSDA in FSDA\_tmp_helpfiles"


echo Dear user, >README.txt
echo since you are reading this file, something went wrong during the installation of FSDA Toolbox helpfiles, probably due to access rights restrictions set on your user account. >>README.txt
echo In this same folder ( _tmp_helpfiles ) you can find a folder named -FSDA- : please move it under the folder -help- of your MATLAB installation path. >>README.txt
echo . >>README.txt
echo Enjoy FSDA. >>README.txt
echo .  >>README.txt

rem echo "dopo readme"
rem pause

rem echo "mathpath " %mathpath%
set matdocroot=%mathpath%\help
rem echo "matdocroot " %matdocroot%


rem move FSDA %matdocroot%\.

mkdir %matdocroot%\FSDA
xcopy FSDA %matdocroot%\FSDA /S /E

cd ..
rem set tmphelp=%where%\FSDA\_tmp_helpfiles\FSDA

rem IF EXIST %tmphelp%  (
rem move non riuscita
rem echo "move non riuscita"
rem ) ELSE (
rem rmdir /S /Q %where%\FSDA\_tmp_helpfiles
rem )


set mathelp=%matdocroot%\FSDA\demos

IF EXIST %mathelp% (
rmdir /S /Q %where%\FSDA\_tmp_helpfiles
) ELSE (
rem echo "move non riuscita"
)

set mathpath=%mathpath:~1,-1%
set where=%where:~1,-1%

"%mathpath%\bin\matlab.exe" -nodesktop -r " builddocsearchdb ('%where%\FSDA\helpfiles\pointersHTML');quit"

rem echo "dopo builddocserachdb"

del /f mgmhlpR7.bat
(goto) 2>nul & del "%~f0"
