#/bin/bash
set -x
pths=`cygpath -a -u ${WORKSPACE}/EXAMPLES_test`
pths_matlab=`cygpath -a -w ${WORKSPACE}/EXAMPLES_test`
wksp_matlab=`cygpath -a -w ${WORKSPACE}`
pth_fsda=`cygpath -a -w ${WORKSPACE}/FSDA`

/usr/bin/find $pths -name "*.m" >flist

# cat flist 

# MATLAB R2012a

rm -f test_runner2012a.m
rm -f execution_log2012a.txt

echo -e "\n" >>test_runner2012a.m

cat flist | while read func_file
do 

    x=`cygpath -w $func_file`

    #    y="clearvars; try; run('$x'); diary('execution_log.txt'); disp('Execution of $x completed successfully'); diary('off'); exit(0); catch error; diary('execution_log.txt'); disp(['Execution of $x FAILED: ' error.message]); diary('off'); exit(0); end;" 
    # NO_EXITS    
    # clearvars 
    # y="try; run('$x'); diary('execution_log2012a.txt'); disp([datestr(clock) ' - Execution of $x completed successfully']); diary('off'); catch error; diary('execution_log.txt'); disp([datestr(clock) ' - Execution of $x FAILED: ' error.message]); diary('off'); end; clear all;" 	
    y="try; run('$x'); diary('execution_log2012a.txt'); disp([datestr(clock) ' - Execution of $x completed successfully']); diary('off'); catch error; diary('execution_log2012a.txt'); disp([datestr(clock) ' - Execution of $x FAILED: ' error.message]); diary('off'); end;" 	
	
    echo $y >>test_runner2012a.m
	  echo -e "\n" >>test_runner2012a.m
done

echo -e "exit(0);\n" >>test_runner2012a.m

# MATLAB R2014b

rm -f test_runner2014b.m
rm -f execution_log2014b.txt

echo -e "\n" >>test_runner2014b.m

cat flist | while read func_file
do 

    x=`cygpath -w $func_file`

    #    y="clearvars; try; run('$x'); diary('execution_log.txt'); disp('Execution of $x completed successfully'); diary('off'); exit(0); catch error; diary('execution_log.txt'); disp(['Execution of $x FAILED: ' error.message]); diary('off'); exit(0); end;" 
    # NO_EXITS    
    # clearvars 
    # y="try; run('$x'); diary('execution_log2014b.txt'); disp([datestr(clock) ' - Execution of $x completed successfully']); diary('off'); catch error; diary('execution_log.txt'); disp([datestr(clock) ' - Execution of $x FAILED: ' error.message]); diary('off'); end; clear all;" 	
    y="try; run('$x'); diary('execution_log2014b.txt'); disp([datestr(clock) ' - Execution of $x completed successfully']); diary('off'); catch error; diary('execution_log2014b.txt'); disp([datestr(clock) ' - Execution of $x FAILED: ' error.message]); diary('off'); end;" 	
	
    echo $y >>test_runner2014b.m
	  echo -e "\n" >>test_runner2014b.m
done

echo -e "exit(0);\n" >>test_runner2014b.m

# MATLAB R2015a

rm -f test_runner2015a.m
rm -f execution_log2015a.txt

echo -e "\n" >>test_runner2015a.m

cat flist | while read func_file
do 

    x=`cygpath -w $func_file`

    #    y="clearvars; try; run('$x'); diary('execution_log.txt'); disp('Execution of $x completed successfully'); diary('off'); exit(0); catch error; diary('execution_log.txt'); disp(['Execution of $x FAILED: ' error.message]); diary('off'); exit(0); end;" 
    # NO_EXITS    
    # clearvars 
    y="try; run('$x'); diary('execution_log2015a.txt'); disp([datestr(clock) ' - Execution of $x completed successfully']); diary('off'); catch error; diary('execution_log2015a.txt'); disp([datestr(clock) ' - Execution of $x FAILED: ' error.message]); diary('off'); end;" 	
	
    echo $y >>test_runner2015a.m
	  echo -e "\n" >>test_runner2015a.m
done

echo -e "exit(0);\n" >>test_runner2015a.m


# '/cygdrive/C/Program Files/MATLAB/R2015a/bin/matlab'
# '/cygdrive/c/Program Files/MATLAB/R2013b/bin/matlab'

if [ $TEST_2012a == "YES" ]; then
'/cygdrive/c/Program Files/MATLAB/R2012a/bin/matlab' -nodisplay -nosplash -noFigureWindows -minimize -wait -r "addpath('$wksp_matlab'); $addpath('$pths_matlab'); addpath(genpath('$pth_fsda')); test_runner2012a"
fi 

if [ $TEST_2014b == "YES" ]; then
'/cygdrive/c/Program Files/MATLAB/R2014b/bin/matlab' -nodisplay -nosplash -noFigureWindows -minimize -wait -r "addpath('$wksp_matlab'); $addpath('$pths_matlab'); addpath(genpath('$pth_fsda')); test_runner2014b"
fi

if [ $TEST_2015a == "YES" ]; then
'/cygdrive/C/Program Files/MATLAB/R2015a/bin/matlab' -nodisplay -nosplash -noFigureWindows -minimize -wait -r "addpath('$wksp_matlab'); $addpath('$pths_matlab'); addpath(genpath('$pth_fsda')); test_runner2015a"
fi