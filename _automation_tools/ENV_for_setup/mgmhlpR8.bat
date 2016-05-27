@echo off

set where=%1

set mathpath=%2

set drv=%where:~1,2%

%drv%

cd "%where%"

rem echo "SONO IN mgmhlpR8 . Il parametro vale " %where% " il drv " %drv%

rem pause

cd FSDA\helpfiles
rename FSDA FSDAtomove
mkdir FSDA
cd FSDAtomove
move helptoc.xml ..\FSDA
move fsda_product_page.html ..\FSDA

cd ..

move FSDAtomove ..\FSDA

cd ..

set matdocroot=%mathpath%\help

move FSDA %matdocroot%\.
pause
%matpath%\bin\matlab.exe -r "builddocsearchdb %matdocroot%\FSDA "
pause

