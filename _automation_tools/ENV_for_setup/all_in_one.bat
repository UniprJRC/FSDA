@echo off
rem set /p SOURCE="Inserisci il path assoluto dell'archivio FSDA --> "
set TOOL=%cd%
cd ..\..
set SOURCE=%cd%

echo "CREAZIONE COPIA TEMPORANEA FSDA ........."
pause
set TEMP="C:\cygwin64\tmp\FSDA_package"
IF EXIST %TEMP% rmdir /S /Q %TEMP%

mkdir %TEMP%

rem IF EXIST %TEMP%\FSDA.OLD rmdir /S /Q %TEMP%\FSDA.old
rem IF EXIST %TEMP%\FSDA move /Y %TEMP%\FSDA %TEMP%\FSDA.old

mkdir %TEMP%\FSDA

Xcopy /E /H /EXCLUDE:%TOOL%\exclude.txt "%SOURCE%"\* %TEMP%\FSDA\.



rem set PWD=%~dp0

cd %TEMP%

rem cd FSDA\helpfiles\FSDA
rem mkdir ..\XML
rem mkdir ..\XML\helpsearch 
rem copy *.xml ..\XML\.
rem copy helpsearch\* ..\XML\helpsearch\.

rem echo "ARCHIVIO jar PER helpfile..............."
rem pause
rem "C:\Program Files\Java\jdk1.7.0_51\bin\jar" cvf ..\XML\help.jar *
rem cd ..
rem rename FSDA FSDAR8
rem rename XML FSDAR7

rem cd ..\..

rem transform_help.exe %PWD%\FSDA\helpfiles\FSDA\ %PWD%\FSDA\helpfiles\FSDAR8\ helptoc.xml %PWD%\template.html

rem rmdir /s /q "%PWD%\FSDA\helpfiles\FSDA"

copy %TOOL%\mgmhlpR7.bat FSDA\.
copy %TOOL%\mgmhlpR8.bat FSDA\.

rem copy brushFAN.mlappinstall FSDA\.
rem copy brushRES.mlappinstall FSDA\.
rem copy brushROB.mlappinstall FSDA\.

echo " mgmhlpRx.bat copiati""
pause

copy %TOOL%\licence.txt .
copy "%TOOL%\before inst.txt" .
copy "%TOOL%\after inst.txt" .
copy %TOOL%\logo.ico .
copy %TOOL%\fsda_black_trasp_300dpi_recol.bmp .
copy %TOOL%\FSDA_logo_trasp_58.bmp .
copy %TOOL%\FSDAscript.iss .
copy %TOOL%\mycode-signing.p12 .


echo "MODIFICA data setup................"
pause

set miadata=%date:/=-%

"C:\cygwin64\bin\sed.exe"  -i "s/AppVerName=FSDA toolbox Release 2017a .*/AppVerName=FSDA toolbox Release 2017a (%miadata%) for MATLAB (not earlier than R2009b) /" FSDAscript.iss

"C:\cygwin64\bin\unix2dos.exe" FSDAscript.iss

echo "COMPILAZIONE setup................"
pause

"C:\Program Files (x86)\Inno Setup 5\Compil32.exe" /cc FSDAscript.iss

echo "----> FSDAtoolbox_for_MATLAB-setup.exe  generato in " %TEMP% " !!!!!!"
pause

"C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\signtool.exe" sign /f .\mycode-signing.p12 /p scacchi2016  %TEMP%\FSDAtoolbox_for_MATLAB-setup.exe

echo "FSDAtoolbox_for_MATLAB-setup.exe SIGNED !"
pause 

echo "Generazione tar package per Linux"
pause

del %TEMP%\FSDA\mgmhlpR7.bat
del %TEMP%\FSDA\mgmhlpR8.bat

copy %TOOL%\setupLINUX.sh .

"C:\cygwin64\bin\tar.exe" cvf FSDA.tar FSDA setupLINUX.sh
"C:\cygwin64\bin\gzip.exe" FSDA.tar
echo "---> FSDA.tar.gz generato in " %TEMP% " !!!!!!!!!!!!!!"
pause


echo "FTP  MR_webserver ........"
pause

ftp -s:%TOOL%\ftp.txt
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
