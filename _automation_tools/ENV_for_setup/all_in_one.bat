@echo off
set /p SOURCE="Inserisci il path assoluto dell'archivio FSDA --> "

if exist "%SOURCE%" goto ok 
if not exist "%SOURCE%" goto err

:ok
echo "RIMOZIONE O RINOMINA COPIA LOCALE FSDA ........."
IF EXIST FSDA.OLD rmdir /S /Q FSDA.old
IF EXIST FSDA move /Y FSDA FSDA.old
mkdir FSDA
echo "COPIA NUOVO ARCHIVIO FSDA ........."
pause
Xcopy /E /H /EXCLUDE:exclude.txt "%SOURCE%"\* FSDA\.

set PWD=%~dp0

cd FSDA\helpfiles\FSDA
mkdir ..\XML
mkdir ..\XML\helpsearch 
copy *.xml ..\XML\.
copy helpsearch\* ..\XML\helpsearch\.

echo "ARCHIVIO jar PER helpfile..............."
pause
"C:\Program Files\Java\jdk1.7.0_51\bin\jar" cvf ..\XML\help.jar *
cd ..
rename FSDA FSDAR8
rename XML FSDAR7
rem mkdir FSDAR8

rem Xcopy /E /H ..\..\new_helps\* FSDAR8\.

cd ..\..

rem transform_help.exe %PWD%\FSDA\helpfiles\FSDA\ %PWD%\FSDA\helpfiles\FSDAR8\ helptoc.xml %PWD%\template.html

rem rmdir /s /q "%PWD%\FSDA\helpfiles\FSDA"

copy mgmhlpR7.bat FSDA\.
copy mgmhlpR8.bat FSDA\.
rem copy brushFAN.mlappinstall FSDA\.
rem copy brushRES.mlappinstall FSDA\.
rem copy brushROB.mlappinstall FSDA\.

echo " mgmhlpRx.bat copiati""
pause

echo "MODIFICA data setup................"
pause

"C:\cygwin64\bin\sed.exe"  -i "s/AppVerName=FSDA toolbox for MATLAB .*/AppVerName=FSDA toolbox for MATLAB Release 3.0 (%date%)/" FSDAscript.iss

"C:\cygwin64\bin\unix2dos.exe" FSDAscript.iss

echo "COMPILAZIONE setup................"
pause

"C:\Program Files (x86)\Inno Setup 5\Compil32.exe" /cc FSDAscript.iss

echo "ECCO FSDAtoolbox_for_MATLAB-setup.exe  !!!!!!"
pause

rem "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\signtool.exe" sign /f "E:\usr_MATLAB\FSDA_setup\CERT\key.pfx" /d "FSDA Toolbox" /du "http://www.riani.it" /t "http://timestamp.verisign.com/scripts/timestamp.dll" "FSDAtoolbox_for_MATLAB-setup.exe"

rem echo "FSDAtoolbox_for_MATLAB-setup.exe SIGNED !"
rem pause 

echo "Generazione tar package per Linux"
pause

del FSDA\mgmhlpR7.bat
del FSDA\mgmhlpR8.bat

"C:\cygwin64\bin\tar.exe" cvf FSDA.tar FSDA setupLINUX.sh
"C:\cygwin64\bin\gzip.exe" FSDA.tar
echo "ECCO FSDA.tar.gz !!!!!!!!!!!!!!"
pause
exit

echo "FTP  MR_webserver ........"
pause

ftp -s:ftp.txt
pause


echo "RICORDA INVIO MAIL A enrico.rossi@ext.ec.europa.eu"
rem echo "FTP  MASTROBUONO server ........"
pause

rem ftp -s:ftpIPSCsite.txt
pause


exit 0

:err
echo "ERRORE !!!!! IL FOLDER " %SOURCE% " NON ESISTE !!"
pause

exit 1
