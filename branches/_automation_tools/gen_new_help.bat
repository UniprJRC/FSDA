@echo off
echo "PROCEDURA per generare nuovi helpfiles REL8 "
echo "<ENTER> per continuare"
pause

set TOOL=%cd%
cd ..
set HELP=%cd%\helpfiles\FSDA
cd %TOOL%\ENV_for_SETUP

.\transform_help.exe %HELP%\ .\ helptoc.xml .\template.html

del /Q indice.js 
del /Q helpfuncbycat.xml
del /Q helpindex.xml
del /Q helptoc.xml

move /Y .\alphabetical.html %HELP%\.
move /Y .\function_reference.html %HELP%\.

echo " RICORDA di fare COMMIT per alphabetical.html e function_reference.html "
echo " <ENTER> per finire "
pause