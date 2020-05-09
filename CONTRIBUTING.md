# Contributing to FSDA

## How do I contribute to FSDA? <a name="toc"></a>


* [Use This Guide](#introduction)
* Ask or Say Something? 
  * [Request Support](#request-support)
  * [Report an Error or Bug](#report-an-error-or-bug)
  * [Request a Feature](#request-a-feature)
* Make Something? 
  * [Project Setup](#project-setup)
  * [Contribute Documentation](#contribute-documentation)
  * [Generate the doc of your own FSDA functions](#generate-the-documentation-of-your-own-FSDA-functions)
  * [Contribute Code](#contribute-code)
* Manage Something 
  * [Label Issues](#label-issues)
  * [Running unit tests](#running-unit-tests)
*  [Join the Project Team](#join-the-project-team)  

## Introduction

Thank you so much for your interest in contributing! All types of contributions are encouraged and valued. See the [table of contents](#toc) for different ways to help and details about how this project handles them!

Please make sure to read the relevant section before making your contribution! It will make it a lot easier for us maintainers to make the most of it and smooth out the experience for all involved. 

The [FSDA Team](#join-the-project-team) looks forward to your contributions. 

## Request Support

If you have a question about this project, how to use it, or just need clarification about something:

* Open an Issue at https://github.com/UniprJRC/FSDA/issues
* Provide as much context as you can about what you're running into.
* Provide platform versions (Windows, Mac, etc), MATLAB version and FSDA version an in general what seems relevant. 

Once it's filed:

* The project team will [label the issue](#label-issues).
* Someone will try to have a response soon.


## Report an Error or Bug

We'd love to accept your patches! 
If you run into an error or bug with the project:

* Open an Issue at Open an Issue at https://github.com/UniprJRC/FSDA/issues
* Include *reproduction steps* that someone else can follow to recreate the bug or error on their own.
* Provide platform versions (Windows, Mac, etc), MATLAB version and FSDA version an in general what seems relevant. 

Once it's filed the project team will [label the issue](#label-issues).

## Request a Feature

If the project doesn't do something you need or want it to do:

* Open an Issue at https://github.com/UniprJRC/FSDA/issues
* Provide as much context as you can about what you're running into.
* Please try and be clear about why existing features and alternatives would not work for you.

Once it's filed the FSDA team will [label the issue](#label-issues).


## Project Setup

So you wanna contribute some code! That's great! This project uses GitHub Pull Requests (PR) to manage contributions, so [read up on how to fork a GitHub project and file a PR](https://guides.github.com/activities/forking) if you've never done it before.

If this seems like a lot or you aren't able to do all this setup, you might also be able to [edit the files directly](https://help.github.com/articles/editing-files-in-another-user-s-repository/) without having to do any of this setup. Yes, [even code](#contribute-code).


## Contribute Documentation

Documentation is a super important, critical part of this project. Docs are how we keep track of what we're doing, how, and why. It's how we stay on the same page about our policies. And it's how we tell others everything they need in order to be able to use this project -- or contribute to it. So thank you in advance.
The FSDA documentation can be found in the supplementary software section of MATLAB help system. A copy of the documentation can be found at the web address below
[http://rosa.unipr.it/FSDA/guide.html](http://rosa.unipr.it/FSDA/guide.html)

![](http://rosa.unipr.it/FSDA/demos/FSDAentrypage.jpg)

Note that this documentation has been automatically created from the .m code files therefore there is a one to one correspondence between what is in the .m file and what is in the associated .html file. See our function [publishFS.m](http://rosa.unipr.it/FSDA/publishFS.html) for further details and section below [Generate the documentation of your own FSDA functions](#generate-the-documentation-of-your-own-FSDA-functions)

Documentation contributions of any size are welcome! Feel free to file a PR even if you're just rewording a sentence to be more clear, or fixing a spelling mistake!

To contribute documentation:

* Edit or add any relevant documentation.
* Make sure your changes are formatted correctly and consistently with the rest of the documentation.
* Re-read what you wrote, and run a spellchecker on it to make sure you didn't miss anything.
* Go to https://github.com/UniprJRC/FSDA/pulls and open a new pull request with your changes.

Before sending your pull requests, make sure you followed the [Code of Conduct](CODE_OF_CONDUCT.md).


Once you've filed the PR:

* One or more maintainers will use GitHub's review feature to review your PR. If your PR gets accepted, it will be marked as such, and merged into the `latest` branch soon after. Your contribution will be distributed to the masses next time the FSDA maintainers [tag a release](https://github.com/UniprJRC/FSDA/releases)

## Generate the documentation of your own FSDA functions

 It is customary for a MATLAB user to document new functions in the head of the .m file. Only rarely the user is prepared to duplicate the effort and work on the corresponding .html documentation file. This is understandable, since the complete integration of new .html files in the standard MATLAB documentation system is not facilitated by built-in tools. In order to help the user in this time consuming but valuable task, FSDA provides some tools, which should be used in the following order:
* [publishFS.m](http://rosa.unipr.it/FSDA/publishFS.html). This is a parser that generates the .html documentation page of a structured .m file. 
The head of the .m file must contain documentation written in compliance with the syntactical rules written in the documentation of publishFS. The syntactical rules are rather intuitive, but are many and not easy to remember if applied episodically. For this reason, it may be a good practice to write documentation starting from an existing FSDA.m file, and refer to the documented syntactical rules just when the parsing errors and warnings seem difficult to interpret. 
The tail of the .m file (i.e. any line after the end function statement) must contain a tag that identifies the category to which the new function is thought to belong. If the tag is not present, the parser does not generate the documentation page. 
* [makecontentsfileFS.m](http://rosa.unipr.it/FSDA/makecontentsfileFS.html). This function generates [personalized .contents files](https://github.com/UniprJRC/FSDA/blob/master/Contents.m) of the functions given in a folder and selected subfolders.  The files to include inside the contents can be filtered according to their filename or content makecontentsFS also returns an output structure containing the list of the files together with information on their location, creation dates, and so on.
The .m function publishFSallFiles calls routine publishFS for the list of files generated by makecontentsFS 
* [publishFunctionAlpha.m](http://rosa.unipr.it/FSDA/publishFunctionCate.html) and [publishFunctionCate.m](http://rosa.unipr.it/FSDA/publishFunctionAlpha.html). These functions generate the [categorical](http://rosa.unipr.it/FSDA/function-cate.html) and [alphabetical](http://rosa.unipr.it/FSDA/function-alpha.html) index pages of the documentation system, starting from the list generated by makecontentsfileFS. An option of  publishFunctionAlpha.m enables us to obain a .txt file (named function-alpha.txt) which contains the names of all files present indexed by makecontentsfileFS separated by commas. In our HTML page automatically created by our parser publishFS we have included a javascript which calls function-alpha.txt and automatically includes a navigation bar to previous and next file in alphabetical order. 
* [publishBibliography.m](http://rosa.unipr.it/FSDA/publishBibliography.html) This functions generates page [bibliography](http://rosa.unipr.it/FSDA/bibliography.html) starting from the output of routine publishFSallFiles. 

## Contribute Code

We like code commits a lot! They're super handy, and they keep the project going and doing the work it needs to do to be useful to others. Code contributions of just about any size are acceptable!

The main difference between code contributions and documentation contributions is that contributing code requires inclusion of relevant tests for the code being added or changed. Contributions without accompanying tests will be held off until a test is added, unless the maintainers consider the specific tests to be either impossible, or way too much of a burden for such a contribution.


If you have improvements to FSDA, send us your pull requests! For those just getting started, Github has a [how to](https://help.github.com/articles/using-pull-requests/).

FSDA team members will be assigned to review your pull requests. Once the
pull requests are approved and pass continuous integration checks, a FSDA
team member will apply `ready to pull` label to your change.  After
the change has been submitted internally, your pull request will be merged
automatically on GitHub.

To contribute code:

* [Set up the project](#project-setup).
* Make any necessary changes to the source code.
* Include any [additional documentation](#contribute-documentation) the changes might need.
* Write tests that verify that your contribution works as expected.
* If your PR gets accepted, it will be marked as such, and merged into the `latest` branch soon after. Your contribution will be distributed in the toolbox next time the maintainers [tag a release](#tag-a-release)



## Label Issues

One of the most important tasks in handling issues is labeling them usefully and accurately. All other tasks involving issues ultimately rely on the issue being classified in such a way that relevant parties looking to do their own tasks can find them quickly and easily.
Below we can find guidelines about how to label the issues.

Label | Apply When | Notes
--- | --- | ---
`bug` | Cases where the code (or documentation) is behaving in a way it wasn't intended to. | If something is happening that surprises the *user* but does not go against the way the code is designed, it should use the `enhancement` label.
`critical` | Added to `bug` issues if the problem described makes the code completely unusable in a common situation. |
`documentation` | Added to issues or pull requests that affect any of the documentation for the project. | Can be combined with other labels, such as `bug` or `enhancement`.
`duplicate` | Added to issues or PRs that refer to the exact same issue as another one that's been previously labeled. | Duplicate issues should be marked and closed right away, with a message referencing the issue it's a duplicate of (with `#123`)
`enhancement` | Added to [feature requests](#request-a-feature), PRs, or documentation issues that are purely additive: the code or docs currently work as expected, but a change is being requested or suggested. |
`help wanted` | Applied by [Committers](#join-the-project-team) to issues and PRs that they would like to get outside help for. Generally, this means it's lower priority for the maintainer team to itself implement, but that the community is encouraged to pick up if they so desire | Never applied on first-pass labeling.
`in-progress` | Applied by [Committers](#join-the-project-team) to PRs that are pending some work before they're ready for review. | The original PR submitter should @mention the team member that applied the label once the PR is complete.
`performance` | This issue or PR is directly related to improving performance. |
`refactor` | Added to issues or PRs that deal with cleaning up or modifying the project for the betterment of it. |
`starter` | Applied by [Committers](#join-the-project-team) to issues that they consider good introductions to the project for people who have not contributed before. These are not necessarily "easy", but rather focused around how much context is necessary in order to understand what needs to be done for this project in particular. | Existing project members are expected to stay away from these unless they increase in priority.
`support` | This issue is either asking a question about how to use the project, clarifying the reason for unexpected behavior, or possibly reporting a `bug` but does not have enough detail yet to determine whether it would count as such. | The label should be switched to `bug` if reliable reproduction steps are provided. Issues primarily with unintended configurations of a user's environment are not considered bugs, even if they cause crashes.
`tests` | This issue or PR either requests or adds primarily tests to the project. | If a PR is pending tests, that will be handled through the PR review process(#review-pull-requests)
`wontfix` | Labelers may apply this label to issues that clearly have nothing at all to do with the project or are otherwise entirely outside of its scope/sphere of influence. [Committers](#join-the-project-team) may apply this label and close an issue or PR if they decide to pass on an otherwise relevant issue. | The issue or PR should be closed as soon as the label is applied, and a clear explanation provided of why the label was used. Contributors are free to contest the labeling, but the decision ultimately falls on committers as to whether to accept something or not.



## Running unit tests

After any pull request from a branch or any modification of a file in the code section, a testing session on Travis, CircleCI and Azure Pipelines is automatically triggered unless the comment section of the commit directive contains string "[skip ci]".


## Join the Project Team

If you arrive at the point of writing new functions and documentation pages compliant with the FSDA philosophy, it means that you have enough energies to take part in the FSDA project. In this case, please check our websites (http://fsda.jrc.ec.europa.eu and http://rosa.unipr.it) for open projects and feel free to contact us at FSDA@unipr.it 
