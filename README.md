![GitHub top language](https://img.shields.io/github/languages/top/UniprJRC/FSDA)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/UniprJRC/FSDA)](https://github.com/UniprJRC/FSDA/releases/latest)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/UniprJRC/FSDA)
[![View FSDA on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72999-fsda)
[![Documentation](https://img.shields.io/badge/HTML_Documentation-Mathworks_style-chocolate.svg)](http://rosa.unipr.it/FSDA/guide.html) 

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FUniprJRC%2FFSDA&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
[![CircleCI](https://circleci.com/gh/UniprJRC/FSDA.svg?style=svg)](https://circleci.com/gh/UniprJRC/FSDA)
[![Build Status](https://dev.azure.com/aldocorbellini0395/FSDA/_apis/build/status/UniprJRC.FSDA%20(1)?branchName=master)](https://dev.azure.com/aldocorbellini0395/FSDA/_build/latest?definitionId=2&branchName=master)
[![MATLAB](https://github.com/UniprJRC/FSDA/actions/workflows/ci.yml/badge.svg)](https://github.com/UniprJRC/FSDA/actions/workflows/ci.yml)

[![codecov](https://codecov.io/gh/UniprJRC/FSDA/branch/master/graph/badge.svg)](https://codecov.io/gh/UniprJRC/FSDA)
[![GitHub contributors](https://img.shields.io/github/contributors/UniprJRC/FSDA)](https://github.com/UniprJRC/FSDA/graphs/contributors)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/UniprJRC/FSDA/graphs/commit-activity)


[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FSDA&project=FSDA.prj)

# [Flexible Robust Statistics Data Analysis](https://github.com/UniprJRC/FSDA/)


## NEWS: new features (June 2023)
1. Now FSDA supports the new [buildtool](https://github.com/UniprJRC/FSDA/tree/master/utilities_help/build) feature (follow the link for more info) to create new releases on GitHub automatically leveraging
the sinergy of the new buildtool functionalities, (see an overview in the [documentation](https://it.mathworks.com/help/matlab/matlab_prog/overview-of-matlab-build-tool.html), available since R2022b)  MATLAB scripts and GitHub Actions.
Releasing a new FSDA release was a manual, multi-step process, that involved a lot of tasks in different environments, now the process runs entirely on GitHub and is consistent and fast.
(Our thanks goes to Jos Martin, Rob Purser, Bensingh Pancras, Andy Campbell, Mark Cafaro et. al. that supported us with the implementation of this new feature)

2. Now FSDA is also availble as a [Docker](https://github.com/UniprJRC/FSDA/tree/master/docker) (follow the link for more info) so if you need to run simulations on a HPC facility on Singularity/Apptainer or you just want to try FSDA with all the features you can follow this link and download a full fledged FSDA docker (yes it works also locally on WSL/WSL2!). 
Once a new release is created, a docker of FSDA is automatically build and can be easily pulled.
(Our thanks go to Jos Martin that helped us a lot on this project).



## FSDA


This project hosts the source code to the [original MATLAB FileExchange project](https://www.mathworks.com/matlabcentral/fileexchange/72999-fsda) and is place of active development.

[![FSDA Toolbox](helpfiles/FSDA/images/githubimgindex.jpg)](http://rosa.unipr.it/FSDA/guide.html)



FSDA Toolbox™ provides statisticians, engineers, scientists, researchers, financial analysts with a comprehensive set of tools to assess and understand their data. Flexible Statistics Data Analysis Toolbox™ software includes functions and interactive tools for analyzing and modeling data, learning and teaching statistics.

The Flexible Statistics Data Analysis Toolbox™ supports a set of routines to develop robust and efficient analysis of complex data sets (multivariate, regression, clustering, ...), ensuring an output unaffected by anomalies or deviations from model assumptions. 

In addition, it offers a rich set interactive graphical tools which enable us to explore the connection in the various features of the different forward plots.

All Flexible Statistics Data Analysis Toolbox™ functions are written in the open MATLAB® language. This means that you can inspect the algorithms, modify the source code, and create your own custom functions.


 For the details about the functions present in FSDA you can browse the categorial and alphabetical list of functions of the toolbox inside MATLAB (once FSDA is installed) or at the web addresses 
 [http://rosa.unipr.it/FSDA/function-cate.html](http://rosa.unipr.it/FSDA/function-cate.html) and [http://rosa.unipr.it/FSDA/function-alpha.html](http://rosa.unipr.it/FSDA/function-alpha.html)




#### FSDA

* Is especially useful in detecting in data potential anomalies (outliers), even when they occur in groups.
Can be used to identify sub-groups in heterogeneous data.
* Extends functionalities in key statistical domains requiring robust analysis (cluster analysis, discriminant analysis, model selection, data transformation).
* Integrates instruments for interactive data visualization and modern exploratory data analysis, designed to simplify the interpretation of the statistical results by the end user.
* Provides statisticians, engineers, scientists, financial analysts a comprehensive set of tools to assess and understand their data.
* Provides practitioners, students and teachers with functions and graphical tools for modeling complex data, learning and teaching statistics.

FSDA is developed for wide applicability. For its capacity to address problems focusing on anomalies in the data, it is expected that it will be used in applications such as anti-fraud, detection of computer network intrusions, e-commerce and credit cards frauds, customer and market segmentation, detection of spurious signals in data acquisition systems, in chemometrics (a wide field covering biochemistry, medicine, biology and chemical engineering), in issues related to the production of official statistics (e.g. imputation and data quality checks), and so on.

For more information see the Wiki page at [https://github.com/UniprJRC/FSDA/wiki](https://github.com/UniprJRC/FSDA/wiki)

#### Ways to familiarize with the FSDA toolbox 


* Run the examples contained in files examples_regression.m or examples_multivariate.m or examples_categorical.m.  Notice that all examples are organized in cells
* Run the GUIs in the FSDA Matlab help pages. 
  For a preview see [http://rosa.unipr.it/FSDA/examples.html](http://rosa.unipr.it/FSDA/examples.html)

  [![FSDA Examples](helpfiles/FSDA/images/githubimgexamples.jpg)](http://rosa.unipr.it/FSDA/examples.html)


* Watch the videos in the Examples section of the FSDA Matlab help pages 
For a preview see [http://rosa.unipr.it/fsda_video.html](http://rosa.unipr.it/fsda_video.html)

* Read section "Introduction to robust statistics" or "Technical introduction to Robust Statistics" in the FSDA Matlab help pages. For a preview see [http://rosa.unipr.it/FSDA/tutorials.html](http://rosa.unipr.it/FSDA/tutorials.html)

  [![FSDA Tutorials](helpfiles/FSDA/images/githubimgtutorials.jpg)](http://rosa.unipr.it/FSDA/tutorials.html)

   


# Installing the toolbox

The installation procedure of FSDA relies on a setup executable (or
Linux bash script) which should execute all the necessary steps
automatically in less than one minute.

MATLAB 64 bit and release from R2012a is required.

The most critical part of the installation concerns the FSDA
documentation system, which consists in a series of HTML files that
follow the typical MATLAB style and are completely integrated inside the
MATLAB documentation system.

In order to achieve this goal from MATLAB release 2012b, it is necessary
that all FSDA html files are located inside the subfolder help of the
main root where MATLAB is installed. For example, if you have release
2017a of MATLAB installed inside C:\\program files\\MATLAB\\R2017a, it
is necessary that our documentation files are inside C:\\program
files\\MATLAB\\R2017a\\help

![](./images/image1.png)

Our SETUP should automatically put our html files inside the appropriate
MATLAB folder. If, due to restrictions on user permissions this was not
possible, you will find the FSDA documentation folder inside subfolder
called \_tmp_helpfiles. If you do not find it, it means that everything
went well and you do not have to do anything, otherwise the subfolder
FSDA in "\_tmp_helpfiles" must be moved manually (possibly by your
system administrator) inside the subfolder help of the main root where
MATLAB is installed (see screenshot above).

Installation notes (details)

In recent years MATLAB undertook several important changes, which made
difficult for us to keep FSDA aligned to the MATLAB system. For example,
from R2012b a Toolstrip replaced menus and toolbars in the MATLAB
Desktop, a gallery of apps was introduced in the desktop and the
documentation system was redesigned, for the first time after years of
stability; with R2014a third party software documentation was moved to a
separate section without filtering and searching possibilities; with
R2014b there was a major update of the MATLAB graphics system that
forced us to revisit most of the FSDA plotting functions; with R2015a
the third party documentation was changed again, and was partly
reintegrated in the MATLAB documentation system. This is why you may
find these installation notes rather intricate: we apologise if this is
going to happen and we invite you to signal problems or bugs that you
think are in contradiction with, or should be reported in, these notes.
With your help, it will be easier for us to make these installation
notes ... superfluous.

1.  FSDA has been fully tested from the release R2012a of MATLAB and
    uses the Statistics toolbox.

2.  FSDA can be installed:

    a.  Automatically with our setup program for Windows platforms. The
        automatic installation updates your MATLAB search path and
        installs the FSDA documentation pages in:

> (docroot folder)\\FSDA, If your release of MATLAB is \<=2012a
>
> (FSDA root folder)\\helpfiles\\FSDA, if your release of MATLAB is
> \>2012a.
>
> To find out where your docroot folder is located it is enough to type
> docroot in the command prompt.

b.  Manually by unpacking the compressed tar file FSDA.tar.gz under a
    folder of your choice (say programs). The search path update can be
    done by running the MATLAB scripts addFSDA2path.m that is located in
    the FSDA root folder. If your release is \>2012a the files present
    inside (FSDA root folder)\\helpfiles\\FSDA must be moved to (docroot
    folder)\\FSDA .

```{=html}
<!-- -->
```
3.  If there are multiple releases of MATLAB installed in your computer,
    our setup program will ask you to **choose to which release the FSDA
    Toolbox has to be associated and where** (under which folder) it has
    to be installed. The search path update and documentation setup are
    modified accordingly. However, if other MATLAB releases are present
    and the user intends to run FSDA also on them, the two steps should
    be completed manually by using the already mentioned addFSDA2path.m
    (see 2b) as follows (now assuming a MS Windows platform installation
    under D:\\programs\\FSDA):

![](./images/image2.png)

4.  If FSDA has been installed properly (in what follows without loss of
    generality we assume, for example, that FSDA has been installed in
    folder D:\\matlab\\FSDA), after the installation the **"Set Path"
    window of MATLAB should include the following FSDA search paths**
    (the last three being introduced with FSDA V3.0 (i.e. starting with
    MATLAB R2015a).

![](./images/image3.png)

5.  If FSDA is installed in MATLAB R2012b or subsequent releases, three
    APPS (brushRES, brushFAN and brushROB) are automatically installed:

![](./images/image4.png)

Remark: if the three APPS have not been automatically installed, you can
easily install them manually by double clicking on the files
brushFAN.mlappinstall, brushRES.mlappinstall and brushROB.mlappinstall
contained in the root folder of FSDA.

![](./images/image5.png)

If FSDA is installed in MATLAB 2012a or earlier the three APPS appear
inside MATLAB Start button\|Toolboxes\|FSDA.

![](./images/image6.png)

> These APPS are graphical user interfaces conceived to demonstrate some
> functionalities of FSDA.

6.  Our setup program, if successfully executed, adds to the "Program
    Files" Windows Menu the entry "FSDA toolbox for MATLAB", including a
    **FSDA uninstall program** that should be used by the user to remove
    an obsolete FSDA release, before an update:

![](./images/image7.png)

7.  Nonetheless, to avoid problems that may occur if FSDA is installed
    with our setup program more than once, the setup program tries to
    locate and remove (with the agreement of the user) previous FSDA
    installations.

8.  If everything went well with an automatic or manual installation,
    when you open MATLAB, the MATLAB "Help" pages should include FSDA
    with all its submenus.

    a.  **The MATLAB "Help" pages should include FSDA**, as shown below:

![](./images/image8.png)

b.  **With the setup installer, three example files named
    "examples_regression.m", "examples_multivariate.m", and
    "examples_categorical.m" should be opened automatically.** These
    files contain a series analysis of several well-known datasets in
    the literature of robust statistics and categorical data analysis
    and have the purpose to let the user familiarize with the toolbox
    (these two files are contained in (main folder of
    FSDA)\\examples)**.**

> ![](./images/image9.png)

c.  **FSDA should appear among the installed "Toolboxes" in the MATLAB
    "Start Menu" (only for MATLAB releases before R2012b)**

![](./images/image10.png)

d.  **For MATLAB R2012b to R2014b installations, the html help files can
    be found in the Supplemental Software tab** which appears at the
    bottom of the Doc Center home page (see screenshot below)**.**

![](./images/image11.png)

**For MATLAB 2015a-2016b installations, the html help files can be found
in the Supplemental Software link**, (see screenshot below)

**Screenshot for MATLAB 2016b**

![](./images/image12.png)

**Screenshot for MATLAB 2017a-2018a**

![](./images/image13.png)

**Screenshot for MATLAB 2019a**

![](./images/image14.png)

e.  **For MATLAB installations earlier than 2012b, the documentation is
    located in the same place as all the other official Mathworks
    toolboxes** (see bottom panel of screenshot below):

![](./images/image15.png)

Independently from MATLAB version you use, once you click on the link
FSDA Toolbox you should reach our main documentation page (see
screenshot below)

![](./images/image8.png)

Remark: you can reach our main documentation page also simply typing
docsearchFS in the command prompt

![](./images/image16.png)

From our main documentation page you can go to the Examples page (see
screenshot below),

![](./images/image17.png)

where you can find GUIs, example codes (see screenshot below),

![](./images/image18.png)

and links to videos containing the analysis of selected examples (see
screenshot below).

![](./images/image19.png)

From any point of our documentation system you can go to the "Tutorials"
page (see screenshot below)

![](./images/image20.png)

where you can find several tutorials about robust statistics and dynamic
statistical visualization, transformations.... (see screenshot below).

![](./images/image21.png)

On the other hand, if from the left menu one clicks on "Functions and
Other References" (see screenshot above), it is possible to get the
categorical list of functions present in the toolbox (see screenshot
below).

![](./images/image22.png)

Of course clicking on the button
![](./images/image23.png) it is possible to browse in alphabetical
order the documentation of the 210 functions present inside the FSDA
toolbox (see screenshot below).

![](./images/image24.png)

By clicking on one of these links (for example on tclust, see screenshot
above) it is possible to reach the HTML documentation of the function in
a perfect new MATLAB documentation style (see screenshot below).

![](./images/image25.png)
These HTML documentation pages have been created automatically by our
routines publishFS. Every HTML documentation contains a series of
**Examples** and Related Examples.

The icon ![](./images/image26.png) at the beginning of the line, indicates
that the associated example has been executed and its output has been
captured inside the HTML file. For example, if you click on the first of
the Related Examples (see screenshot below),

![](./images/image27.png)

it is possible to see both the code (note that the code is displayed
inside HTML using typical Matlab colouring) and the output which was
generated (see the two screenshots below).

![](./images/image28.png)

![](./images/image29.png)

In the More About section of our HTML files (see screenshot below), it
is possible find the theoretical background which accompanies a
particular function.

For example, the screenshot below shows what you get in the case of
function tclust.

![](./images/image30.png)

Remark: there is a one to one correspondence between the documentation
contained inside the .m file and the corresponding .html file.

The documentation inside the .m file can be easily accessed from the
command prompt typing help and the name of the function.

For example, the screenshot below shows what you get if you type in the
prompt "help MixSim".

![](./images/image31.png)

Sometimes inside the .m file (especially in the section "More About") we
have added a number of formulae in latex language (see screenshot
below).

![](./images/image32.png)

Clearly all these latex formulae will show up correctly (thanks to
MathJax technology) in the corresponding HTML help page. For example, in
the case of MixSim function, in the command window, by clicking on the
link "Link to the help function",

![](./images/image33.png)

one is redirected to the corresponding HTML documentation page. Here, in
the "More About" section it is possible to see the code in proper
mathematical style.

![](./images/image34.png)

Finally, it is worthwhile to remark that it is possible to go directly
to the HTML documentation page simply typing docsearchFS and the name of
the requested function. For example, in the case of tclust to reach file
tclust.html it is possible to type:

![](./images/image35.png)

Generally, the output of our functions is a structure, which contains
several fields, documented in detail inside the initial part of the .m
function. For example, in the case of tclust inside tclust.m it is
possible to navigate to section Output(see screenshot below):

![](./images/image36.png)

In the corresponding HTML file our parser publishFS.m puts all the
fields of input and output structure inside a HTML table (see screenshot
below):

![](./images/image37.png)

Every subfolder of FSDA contains file contents.m (automatically created
by our routine makecontentsfileFS.m) which contains a series of detailed
information about all the .m files of the folder, which have the
corresponding HTML documentation. For example, the screenshot referred
to the left part of file contents.m inside subfolder "utilities" is
given below.

![](./images/image38.png)

Similarly, inside the main root of FSDA file contents.m lists in
alphabetical order all files present in all subfolders of FSDA, which
have the corresponding HTML page (see screenshot below):

![](./images/image39.png)

## FSDA html documentation files and MATLAB search engine

Particular attention has been devoted by the FSDA team to have all our
HTML files indexed by the old and new search engine of MATLAB. Below we
describe what you should get depending on the MATLAB version you have.

### MATLAB 2015a-2018b

If your version of MATLAB is in the range 2015a-2018b, typing inside the
engine a word you get also the results from the FSDA toolbox. For
example, typing tclust you should automatically have the search
suggestion from the drop down menu which automatically appears (see
screenshot below)

![](./images/image40.png)

and you should be brought directly to the tclust documentation page.

If, for example, you type "concentration step" and you do Refine by
Product and select the FSDA toolbox these are the three instances you
should get.

![](./images/image41.png)

If you put your mouse on the word restreigen you can see from the status
bar that the function is located inside (main FSDA
folder)/helpfiles/pointersHTML, (screenshot of status bar is given
below):

![](./images/image42.png)
Once you click on restreigen you can go to page restreigen.html (see
screenshot below) which is located inside docrootFS/FSDA.

![](./images/image43.png)

From the toolstrip on top you will notice that two instances of
restreigen have been opened.
![](./images/image44.png)
The first on the left is the page which has been indexed by MATLAB
search engine which is located inside

(main FSDA folder)/helpfiles/pointersHTML (see screenshot below):

![](./images/image45.png)
All these syllabus pages have been automatically created by our routine
createFSDAhelp.m. It was necessary to have the intermediate pages
because MATLAB forces these pages to be opened on the iframe on the
right. All these syllabus pages contain a redirect to our final HTML
pages, which are contained inside docroot/FSDA. Files inside
docroot/FSDA are not forced to be opened on the iframe on the right.

### MATLAB 2012b-2014b

If your version of MATLAB is between 2012b-2014b, it is necessary to use
the old MATLAB search engine inside supplemental software (see
screenshot below),

![](./images/image46.png)

also, in this case, passing through the syllabus page contained in (FSDA
root)\\helpfiles\\pointersHTML.

![](./images/image47.png)

It is possible to automatically reach our page LXS.html contained inside
docroot/FSDA.

![](./images/image48.png)

### MATLAB \<=2012a

In MATLAB older or equal than 2012a, there was no distinction between
MATLAB toolboxes and third parties toolboxes (as concerns the
documentation), therefore it is possible to search directly from the
unique official MATLAB engine. For example searching for LXS,

![](./images/image49.png)
the output of the search is again the syllabus page which automatically
redirects to the true HTML page (both are shown in the screenshot
below):

![](./images/image50.png)

Remark 1: note that the old MATLAB browser enables us to see correctly
just 95% of the javascripts which characterize the new MATLAB help
style.

Remark 2: given that in MATLAB versions earlier than 2012b the new
engine lucene did not exist, all the "true HTML files" during the
installation are not moved to folder docroot/FSDA but are left inside
(main root of FSDA)/helpfiles/FSDA.

Remark 3: If you are using a release lower than R2012b and you think
that the MATLAB Help Browser is not producing proper search results for
FSDA functions, check first that in the MATLAB Help Preferences FSDA is
selected, as shown in Figure below:

> ![](./images/image51.JPG)

Remark 4: **From MATLAB R2015a,** **when you search for a given third
party function** **for the first time**, the search results window will
display a yellow message warning that a toolbox in the path does not
have the proper documentation index file. The window and message
produced when attempting to search for documentation about \`FSR'
function, are shown here:

![](./images/image52.png)

> Only after clicking on the Build Index button you should start getting
> the desired documentation. You will be informed of the successful
> update of the search database with this message:
>
> ![](./images/image53.png)
>
> **If,** instead of this message and instead of receiving back the
> desired results, **you receive an error message** such as this one

![](./images/image54.png)

> again it is likely that you have installed FSDA in a location without
> proper permissions and, thus, the index building operation (i.e. the
> builddocsearchdb command) could not update the search database. The
> only solution in this case is to obtain the writing permissions or to
> change location for FSDA. What you should get if the search is
> successful, is something like the following:

![](./images/image55.png)

Therefore, if by some obscure reason you cannot find our HTML files
using old or new (lucene) engine, it might be necessary to run again
buildocsearchdb. For example assuming that FSDA main folder is
D:\\PACKAGES\\FSDA, then it is necessary to run:

![](./images/image56.png)
If you think that something not described in these
notes went wrong please do not hesitate to send an e-mail to

<FSDA@unipr.it>


