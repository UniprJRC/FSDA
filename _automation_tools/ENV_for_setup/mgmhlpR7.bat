@echo off

set where=%1

set drv=%where:~1,2%

%drv%

cd "%where%"


rem echo "SONO IN mgmhlpR7 . Il parametro vale " %where%

cd FSDA\helpfiles

rem rmdir /S /Q FSDAR8

rem rename FSDAR7 FSDA