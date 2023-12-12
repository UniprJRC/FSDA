![](RackMultipart20231107-1-an42vi_html_5b83af8e3b949f63.png)

# Installation notes (highlights)

The installation procedure of FSDA relies on a setup executable (or Linux bash script) which should execute all the necessary steps automatically in less than one minute.

MATLAB 64 bit and release from R2012a is required.

The most critical part of the installation concerns the FSDA documentation system, which consists in a series of HTML files that follow the typical MATLAB style and are completely integrated inside the MATLAB documentation system.

In order to achieve this goal from MATLAB release 2012b, it is necessary that all FSDA html files are located inside the subfolder help of the main root where MATLAB is installed. For example, if you have release 2017a of MATLAB installed inside C:\program files\MATLAB\R2017a, it is necessary that our documentation files are inside C:\program files\MATLAB\R2017a\help

![](RackMultipart20231107-1-an42vi_html_c57ee9cb3434d38a.png)

Our SETUP should automatically put our html files inside the appropriate MATLAB folder. If, due to restrictions on user permissions this was not possible, you will find the FSDA documentation folder inside subfolder called \_tmp\_helpfiles. If you do not find it, it means that everything went well and you do not have to do anything, otherwise the subfolder FSDA in "\_tmp\_helpfiles" must be moved manually (possibly by your system administrator) inside the subfolder help of the main root where MATLAB is installed (see screenshot above).

# Installation notes (details)

In recent years MATLAB undertook several important changes, which made difficult for us to keep FSDA aligned to the MATLAB system. For example, from R2012b a Toolstrip replaced menus and toolbars in the MATLAB Desktop, a gallery of apps was introduced in the desktop and the documentation system was redesigned, for the first time after years of stability; with R2014a third party software documentation was moved to a separate section without filtering and searching possibilities; with R2014b there was a major update of the MATLAB graphics system that forced us to revisit most of the FSDA plotting functions; with R2015a the third party documentation was changed again, and was partly reintegrated in the MATLAB documentation system. This is why you may find these installation notes rather intricate: we apologise if this is going to happen and we invite you to signal problems or bugs that you think are in contradiction with, or should be reported in, these notes. With your help, it will be easier for us to make these installation notes … superfluous.

1. FSDA has been fully tested from the release R2012a of MATLAB and uses the Statistics toolbox.
2. FSDA can be installed:
  1. Automatically with our setup program for Windows platforms. The automatic installation updates your MATLAB search path and installs the FSDA documentation pages in:

(docroot folder)\FSDA, If your release of MATLAB is \<=2012a

(FSDA root folder)\helpfiles\FSDA, if your release of MATLAB is \>2012a.

To find out where your docroot folder is located it is enough to type docroot in the command prompt.

  1. Manually by unpacking the compressed tar file FSDA.tar.gz under a folder of your choice (say programs). The search path update can be done by running the MATLAB scripts addFSDA2path.m that is located in the FSDA root folder. If your release is \>2012a the files present inside (FSDA root folder)\helpfiles\FSDA must be moved to (docroot folder)\FSDA .
1. If there are multiple releases of MATLAB installed in your computer, our setup program will ask you to **choose to which release the FSDA Toolbox has to be associated and where** (under which folder) it has to be installed. The search path update and documentation setup are modified accordingly. However, if other MATLAB releases are present and the user intends to run FSDA also on them, the two steps should be completed manually by using the already mentioned addFSDA2path.m (see 2b) as follows (now assuming a MS Windows platform installation under D:\programs\FSDA):

![](RackMultipart20231107-1-an42vi_html_9a738e286415ef98.png)

1. If FSDA has been installed properly (in what follows without loss of generality we assume, for example, that FSDA has been installed in folder D:\matlab\FSDA), after the installation the **"Set Path" window of MATLAB should include the following FSDA search paths** (the last three being introduced with FSDA V3.0 (i.e. starting with MATLAB R2015a).

![](RackMultipart20231107-1-an42vi_html_2495dd56d0a72f40.png)

1. If FSDA is installed in MATLAB R2012b or subsequent releases, three APPS (brushRES, brushFAN and brushROB) are automatically installed:

![](RackMultipart20231107-1-an42vi_html_87576d778da020e0.png)

Remark: if the three APPS have not been automatically installed, you can easily install them manually by double clicking on the files brushFAN.mlappinstall, brushRES.mlappinstall and brushROB.mlappinstall contained in the root folder of FSDA.

![](RackMultipart20231107-1-an42vi_html_866ae2c396fa817.png)

If FSDA is installed in MATLAB 2012a or earlier the three APPS appear inside MATLABStart button|Toolboxes|FSDA.

![](RackMultipart20231107-1-an42vi_html_d1ce69345d75a175.png)

These APPS are graphical user interfaces conceived to demonstrate some functionalities of FSDA.

1. Our setup program, if successfully executed, adds to the "Program Files" Windows Menu the entry "FSDA toolbox for MATLAB", includinga **FSDA uninstall program** that should be used by the user to remove an obsolete FSDA release, before an update:

![Shape1](RackMultipart20231107-1-an42vi_html_ede67c66c065e128.gif)

1. Nonetheless, to avoid problems that may occur if FSDA is installed with our setup program more than once, the setup program tries to locate and remove (with the agreement of the user) previous FSDA installations.
2. If everything went well with an automatic or manual installation, when you open MATLAB, the MATLAB "Help" pages should include FSDA with all its submenus.
  1. **The MATLAB "Help" pages should include FSDA** , as shown below:

![](RackMultipart20231107-1-an42vi_html_c8efcc58a7c22e91.png)

  1. **With the setup installer, three example files named "examples\_regression.m", "examples\_multivariate.m", and "**** examples\_categorical.m" **** should be opened automatically.**These files contain a series analysis of several well-known datasets in the literature of robust statistics and categorical data analysis and have the purpose to let the user familiarize with the toolbox (these two files are contained in (main folder of FSDA)\examples)**.**

![](RackMultipart20231107-1-an42vi_html_48944aad1ec214c8.png)

  1. **FSDA should appear among the installed "Toolboxes" in the MATLAB "Start Menu" (only for MATLAB releases before R2012b)**

![](RackMultipart20231107-1-an42vi_html_e2b77f3343b524e.png)

  1. **For MATLAB R2012b to R2014b installations, the html help files can be found in the Supplemental Software tab** which appears at the bottom of the Doc Center home page(see screenshot below) **.**

![](RackMultipart20231107-1-an42vi_html_b38498ac9152431c.png)

**For MATLAB 2015a-2016b installations, the html help files can be found in the Supplemental Software link** , (see screenshot below)

**Screenshot for MATLAB 2016b**

![](RackMultipart20231107-1-an42vi_html_dcdba06f2a9e4abd.png)

**Screenshot for MATLAB 2017a-2018a**

![](RackMultipart20231107-1-an42vi_html_ed52407d37fbbe1e.png)

**Screenshot for MATLAB 2019a**

![](RackMultipart20231107-1-an42vi_html_3650146bb377d430.png)

  1. **For MATLAB installations earlier than 2012b, the documentation is located in the same place as all the other official Mathworks toolboxes** (see bottom panel of screenshot below):

![](RackMultipart20231107-1-an42vi_html_e9f84b956d0e4b87.png)

Independently from MATLAB version you use, once you click on the link FSDA Toolbox you should reach our main documentation page (see screenshot below)

![](RackMultipart20231107-1-an42vi_html_c8efcc58a7c22e91.png)

Remark: you can reach our main documentation page also simply typing docsearchFS in the command prompt

![](RackMultipart20231107-1-an42vi_html_c982c661dc2f70c1.png)

From our main documentation page you can go to the Examples page (see screenshot below),

![](RackMultipart20231107-1-an42vi_html_dbeaa4b28c141589.png)

where you can find GUIs, example codes (see screenshot below),

![](RackMultipart20231107-1-an42vi_html_68d88c4d8022cc2b.png)

and links to videos containing the analysis of selected examples (see screenshot below).

![](RackMultipart20231107-1-an42vi_html_cc9dd45e01b8c1bb.png)

From any point of our documentation system you can go to the "Tutorials" page (see screenshot below)

![](RackMultipart20231107-1-an42vi_html_1900165d778e2b96.png)

where you can find several tutorials about robust statistics and dynamic statistical visualization, transformations…. (see screenshot below).

![](RackMultipart20231107-1-an42vi_html_dc1924fd8b0017cd.png)

On the other hand, if from the left menu one clicks on "Functions and Other References" (see screenshot above), it is possible to get the categorical list of functions present in the toolbox (see screenshot below).

![](RackMultipart20231107-1-an42vi_html_9ca8e7b76f920454.png)

Of course clicking on the button ![](RackMultipart20231107-1-an42vi_html_db1bbc8a09854af4.png) it is possible to browse in alphabetical order the documentation of the 210 functions present inside the FSDA toolbox (see screenshot below).

![](RackMultipart20231107-1-an42vi_html_33f0cd16ba746c16.png)

By clicking on one of these links (for example on tclust, see screenshot above) it is possible to reach the HTML documentation of the function in a perfect new MATLAB documentation style (see screenshot below).

![](RackMultipart20231107-1-an42vi_html_2ac702207bd31ae2.png)

These HTML documentation pages have been created automatically by our routines publishFS. Every HTML documentation contains a series of **Examples** and Related Examples.

The icon ![](RackMultipart20231107-1-an42vi_html_3f684fbd35178e1a.png) at the beginning of the line, indicates that the associated example has been executed and its output has been captured inside the HTML file. For example, if you click on the first of the Related Examples (see screenshot below),

![](RackMultipart20231107-1-an42vi_html_bed0ab9897628303.png)

it is possible to see both the code (note that the code is displayed inside HTML using typical Matlab colouring) and the output which was generated (see the two screenshots below).

![](RackMultipart20231107-1-an42vi_html_1ce3171fd541f743.png)

![](RackMultipart20231107-1-an42vi_html_8b2d10c5129b3a9.png)

In the More About section of our HTML files (see screenshot below), it is possible find the theoretical background which accompanies a particular function.

For example, the screenshot below shows what you get in the case of function tclust.

![](RackMultipart20231107-1-an42vi_html_76e2e332a60cf875.png)

Remark: there is a one to one correspondence between the documentation contained inside the .m file and the corresponding .html file.

The documentation inside the .m file can be easily accessed from the command prompt typing help and the name of the function.

For example, the screenshot below shows what you get if you type in the prompt "help MixSim".

![](RackMultipart20231107-1-an42vi_html_69d9d105956c3440.png)

Sometimes inside the .m file (especially in the section "More About") we have added a number of formulae in latex language (see screenshot below).

![](RackMultipart20231107-1-an42vi_html_ad14bba09924988f.png)

Clearly all these latex formulae will show up correctly (thanks to MathJax technology) in the corresponding HTML help page. For example, in the case of MixSim function, in the command window, by clicking on the link "Link to the help function",

![](RackMultipart20231107-1-an42vi_html_975a781a4af6d182.png)

one is redirected to the corresponding HTML documentation page. Here, in the "More About" section it is possible to see the code in proper mathematical style.

![](RackMultipart20231107-1-an42vi_html_a878e436b02346b5.png)

Finally, it is worthwhile to remark that it is possible to go directly to the HTML documentation page simply typing docsearchFS and the name of the requested function. For example, in the case of tclust to reach file tclust.html it is possible to type:

![](RackMultipart20231107-1-an42vi_html_5945ab9ae49b396f.png)

Generally, the output of our functions is a structure, which contains several fields, documented in detail inside the initial part of the .m function. For example, in the case of tclust inside tclust.m it is possible to navigate to section Output(see screenshot below):

![](RackMultipart20231107-1-an42vi_html_2d4d20385781cb83.png)

In the corresponding HTML file our parser publishFS.m puts all the fields of input and output structure inside a HTML table (see screenshot below):

![](RackMultipart20231107-1-an42vi_html_cea85f86330dd2f7.png)

Every subfolder of FSDA contains file contents.m (automatically created by our routine makecontentsfileFS.m) which contains a series of detailed information about all the .m files of the folder, which have the corresponding HTML documentation. For example, the screenshot referred to the left part of file contents.m inside subfolder "utilities" is given below.

![](RackMultipart20231107-1-an42vi_html_9bd156fe027b1de1.png)

Similarly, inside the main root of FSDA file contents.m lists in alphabetical order all files present in all subfolders of FSDA, which have the corresponding HTML page (see screenshot below):

![](RackMultipart20231107-1-an42vi_html_dd48eff41b323e2c.png)

## FSDA html documentation files and MATLAB search engine

Particular attention has been devoted by the FSDA team to have all our HTML files indexed by the old and new search engine of MATLAB. Below we describe what you should get depending on the MATLAB version you have.

### MATLAB 2015a-2018b

If your version of MATLAB is in the range 2015a-2018b, typing inside the engine a word you get also the results from the FSDA toolbox. For example, typing tclust you should automatically have the search suggestion from the drop down menu which automatically appears (see screenshot below)

![](RackMultipart20231107-1-an42vi_html_6b88da109dbaf944.png)

and you should be brought directly to the tclust documentation page.

If, for example, you type "concentration step" and you do Refine by Product and select the FSDA toolbox these are the three instances you should get.

![](RackMultipart20231107-1-an42vi_html_baaef41d61fb1580.png)

If you put your mouse on the word restreigen you can see from the status bar that the function is located inside (main FSDA folder)/helpfiles/pointersHTML, (screenshot of status bar is given below):

![](RackMultipart20231107-1-an42vi_html_b8f8aec5d27e2d92.png)

Once you click on restreigen you can go to page restreigen.html (see screenshot below) which is located inside docrootFS/FSDA.

![](RackMultipart20231107-1-an42vi_html_444f043bc10bb860.png)

From the toolstrip on top you will notice that two instances of restreigen have been opened. ![](RackMultipart20231107-1-an42vi_html_7609144513964f90.png)

The first on the left is the page which has been indexed by MATLAB search engine which is located inside

(main FSDA folder)/helpfiles/pointersHTML (see screenshot below):

![](RackMultipart20231107-1-an42vi_html_a1c1635c0d5d9bfa.png)

All these syllabus pages have been automatically created by our routine createFSDAhelp.m. It was necessary to have the intermediate pages because MATLAB forces these pages to be opened on the iframe on the right. All these syllabus pages contain a redirect to our final HTML pages, which are contained inside docroot/FSDA. Files inside docroot/FSDA are not forced to be opened on the iframe on the right.

### MATLAB 2012b-2014b

If your version of MATLAB is between 2012b-2014b, it is necessary to use the old MATLAB search engine inside supplemental software (see screenshot below),

![](RackMultipart20231107-1-an42vi_html_7afef5f6eae1fac5.png)

also, in this case, passing through the syllabus page contained in (FSDA root)\helpfiles\pointersHTML.

![](RackMultipart20231107-1-an42vi_html_6becdb0b01e6759a.png)

It is possible to automatically reach our page LXS.html contained inside docroot/FSDA.

![](RackMultipart20231107-1-an42vi_html_f395f6b607674110.png)

### MATLAB \<=2012a

In MATLAB older or equal than 2012a, there was no distinction between MATLAB toolboxes and third parties toolboxes (as concerns the documentation), therefore it is possible to search directly from the unique official MATLAB engine. For example searching for LXS,

![](RackMultipart20231107-1-an42vi_html_fc2e14949e61d482.png)

the output of the search is again the syllabus page which automatically redirects to the true HTML page (both are shown in the screenshot below):

![](RackMultipart20231107-1-an42vi_html_2343029b8fd7eafe.png)

Remark 1: note that the old MATLAB browser enables us to see correctly just 95% of the javascripts which characterize the new MATLAB help style.

Remark 2: given that in MATLAB versions earlier than 2012b the new engine lucene did not exist, all the "true HTML files" during the installation are not moved to folder docroot/FSDA but are left inside (main root of FSDA)/helpfiles/FSDA.

Remark 3: If you are using a release lower than R2012b and you think that the MATLAB Help Browser is not producing proper search results for FSDA functions, check first that in the MATLAB Help Preferences FSDA is selected, as shown in Figure below:

![](RackMultipart20231107-1-an42vi_html_4379e3a2543ad438.jpg)

Remark 4: **From MATLAB R2015a,**** when you search for a given third party function ****for the first time** , the search results window will display a yellow message warning that a toolbox in the path does not have the proper documentation index file. The window and message produced when attempting to search for documentation about `FSR' function, are shown here:

![](RackMultipart20231107-1-an42vi_html_3dd6c7dd2678ec63.png)

Only after clicking on the Build Index buttonyou should start getting the desired documentation. You will be informed of the successful update of the search database with this message:

![](RackMultipart20231107-1-an42vi_html_2c2c799c2a3d36a3.png)

**If,** instead of this message and instead of receiving back the desired results, **you receive an error message** such as this one

![](RackMultipart20231107-1-an42vi_html_6f881732e3296a14.png)

again it is likely that you have installed FSDA in a location without proper permissions and, thus, the index building operation (i.e. the builddocsearchdb command) could not update the search database. The only solution in this case is to obtain the writing permissions or to change location for FSDA. What you should get if the search is successful, is something like the following:

![](RackMultipart20231107-1-an42vi_html_8babc66dd7181820.png)

Therefore, if by some obscure reason you cannot find our HTML files using old or new (lucene) engine, it might be necessary to run again buildocsearchdb. For example assuming that FSDA main folder is D:\PACKAGES\FSDA, then it is necessary to run:

![](RackMultipart20231107-1-an42vi_html_7e03714b9d3e8052.png)If you think that something not described in these notes went wrong

please do not hesitate to send an e-mail to

[FSDA@unipr.it](mailto:FSDA@unipr.it)

Page 5/5