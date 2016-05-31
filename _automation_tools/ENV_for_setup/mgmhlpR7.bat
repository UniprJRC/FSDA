@echo off

set where=%1

set mathpath=%2

set drv=%where:~1,2%

%drv%

cd "%where%"

set where=%where:~1,-1%

rem echo "SONO IN mgmhlpR7 . Il parametro vale " %where%

%mathpath%\bin\matlab.exe -nodesktop -r " builddocsearchdb ('%where%\FSDA\helpfiles\FSDA') "  

rem rmdir /S /Q FSDAR8

rem rename FSDAR7 FSDA