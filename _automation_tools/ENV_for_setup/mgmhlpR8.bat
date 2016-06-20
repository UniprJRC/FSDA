@echo off
set where=%1

set mathpath=%2

set drv=%where:~1,2%

%drv%

cd "%where%"

rem echo "SONO IN mgmhlpR8 . Il parametro vale " %where% " il drv " %drv%

rem %mathpath%\bin\matlab.exe -wait -automation -nodesktop -r " builddocsearchdb ('%where%\FSDA\helpfiles\FSDA') ; quit "

cd FSDA\helpfiles
rename FSDA FSDAtomove
mkdir FSDA
cd FSDAtomove
move helptoc.xml ..\FSDA
move fsda_product_page.html ..\FSDA

cd ..

move FSDAtomove ..\FSDA

cd ..

rem echo "mathpath " %mathpath%
set matdocroot=%mathpath%\help
rem echo "matdocroot " %matdocroot%
rem pause

move FSDA %matdocroot%\.

set mathpath=%mathpath:~1,-1%
set matdocroot=%mathpath%\help

rem pause
rem "%mathpath%"\bin\matlab.exe -nodesktop -r " addpath '%matdocroot%\FSDA' ; builddocsearchdb ('%matdocroot%\FSDA') "
rem pause
del /f mgmhlpR7.bat
(goto) 2>nul & del "%~f0"
