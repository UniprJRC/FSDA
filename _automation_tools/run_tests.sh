#/bin/bash
set -x
pths=`cygpath -a -u ${WORKSPACE}/EXAMPLES_test`
pths_matlab=`cygpath -a -w ${WORKSPACE}/EXAMPLES_test`
wksp_matlab=`cygpath -a -w ${WORKSPACE}`
pth_fsda=`cygpath -a -w ${WORKSPACE}/FSDA`

/usr/bin/find $pths -name "*.m" >flist

# cat flist 

rm -f test_runner.m
rm -f execution_log.txt

echo -e "\n" >>test_runner.m

cat flist | while read func_file
do 

    x=`cygpath -w $func_file`

    #    y="clearvars; try; run('$x'); diary('execution_log.txt'); disp('Execution of $x completed successfully'); diary('off'); exit(0); catch error; diary('execution_log.txt'); disp(['Execution of $x FAILED: ' error.message]); diary('off'); exit(0); end;" 
    # NO_EXITS    
    # clearvars 
    y="try; run('$x'); diary('execution_log.txt'); disp([datestr(clock) ' - Execution of $x completed successfully']); diary('off'); catch error; diary('execution_log.txt'); disp([datestr(clock) ' - Execution of $x FAILED: ' error.message]); diary('off'); end; clear all;" 	
	
    echo $y >>test_runner.m
	echo -e "\n" >>test_runner.m
done

echo -e "exit(0);\n" >>test_runner.m

'/cygdrive/c/Program Files/MATLAB/R2013b/bin/matlab' -nodisplay -nosplash -noFigureWindows -minimize -wait -r "addpath('$wksp_matlab'); $addpath('$pths_matlab'); addpath(genpath('$pth_fsda')); test_runner"
