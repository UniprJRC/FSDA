@echo off

set where=%1

set matexe=%2

set drv=%where:~1,2%

%drv%

cd "%where%"

rem echo "SONO IN mgmhlpR8 . Il parametro vale " %where% " il drv " %drv%

rem pause

cd FSDA\helpfiles

rmdir /S /Q FSDAR7

rename FSDAR8 FSDA

set matpat=%matexe:bin\matlab.exe=%

set matpatN=%matpat:\=\\%

set matpatF=%matpatN:"=%

set final=help\\documentation-center.html

rem echo var home = 'file:///%matpatF%%final%';  >FSDA\home.js
echo var home = 'fsda_product_page.html' ; >FSDA\home.js