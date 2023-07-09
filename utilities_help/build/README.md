# Introduction
This folder contains the various build action scripts needed to create the 
MLTBX toolbox file for the FSDA toolbox. These functions can be called on 
both a local machine and in a github action to build the toolbox.

They are used by the `buildfile.m` function in the root of the repository
so that invokations of `buildtool doc` and `buildtool toolbox` generate doc 
searches and the toolbox build respectively.

# Github Release Procedure
To trigger a completely new release of the toolbox on `github` using the 
`upload-artifact.yml` workflow from the main branch of code do the following

1. Clone a completely clean version of the codebase from github using 
```
    git clone https://github.com/UniprJRC/FSDA.git
```
2. In that clean repo create the tag for the new release of the FSDA toolbox
```
    git tag -a -m "Release comment" 8.7.1.0
```
3. Push the new tags to github
```
    git push origin 8.7.1.0
```

####Update (June 2023)

A new possibility is to use `buildtool` feature inside MATLAB: to create a new release type the following in the command window
```
>> buildtool releaseToGithub(Version="8.7.1.0", Comment="This is a comment for the release")
```



This push of a git tag will trigger the `upload-artifact` workflow, since it specifies:
``` yaml 
on:
  push:
    # Sequence of patterns matched against refs/tags
    tags:
    - '*' # Push events to matching *, i.e. v1.0 or 18.6.1
```

This will trigger 3 different dependent workflows (`build-docsearch-v3`, `archive-build-artifacts`, 
`build-docker-container`) that run in that order. The first workflow simply builds pre-22a doc search folders for
the toolbox. The second is the toolbox packaging workflow which download the doc search from the first workflow,
builds the toolbox using the `buildtool` via `run-build`, and the final builds a docker container as described in `docker\README.md`.

These workflows can be seen in the github actions area of the repository, and if successful a new release
will be available on the [repo release page](https://github.com/UniprJRC/FSDA/releases). There will be the usual 
source code assets along with the `FSDA.mltbx` file for customer installation.