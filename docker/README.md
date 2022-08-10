# Introduction
This folder builds a docker container with all the pre-requisites for the [FSDA Toolbox](http://rosa.unipr.it/fsda.html). In addition there are automated github actions to build the same docker container and store it in  `ghcr.io` so that users can simply pull the pre-built container.

## Getting Pre-Built Images
An example of how to `pull` a specific instance (both FSDA version and MATLAB version) is given below:

```
docker pull ghcr.io/uniprjrc/fsda:8.5.29-r2022a
```

Fo a full list of the available tags for these docker containers please see the [repository packages](https://github.com/UniprJRC/FSDA/pkgs/container/fsda).

## Running Container
This container is built directly [`FROM mathworks/matlab:latest`](https://hub.docker.com/r/mathworks/matlab) and has all the same command-line arguments as listed on its help page.

## Build from Repo
To build this container clone the repository and execute

```
docker build -t fsda --build-arg FSDA_RELEASE=${REQUIRED_RELEASE} ./docker
```

You can specify previous versions of FSDA by modifying the `FSDA_RELEASE=X.Y.Z` build argument (see [available releases](https://github.com/UniprJRC/FSDA/releases)).

By default this will build against the `latest` MATLAB docker container but you can specify a particular release by setting the `build-arg` `MATLAB_DOCKER_TAG` to something like `r2022a` (or [any other tag](https://hub.docker.com/r/mathworks/matlab) supported for the docker container).

If you wish to install from your own forked repo of [`UniprJRC/FSDA`](https://github.com/UniprJRC/FSDA) then please set the `build-arg` `SRC_REPO` to the name of your own github repository.

# GitHub Automated Build
The automated build procedure is driven when a specific tag of the form `X.Y.Z` is pushed to the repository. In addition you can trigger a manual build of a particular release using the GitHub CLI (or browse to the action on github). The CLI trigger action is something like

```
gh workflow run -F fsda-release=8.5.29 -F matlab-release=r2022a build-docker-container
```

This will trigger a build that results in a container called `ghcr.io/uniprjrc/fsda:8.5.29-r2022a`.
