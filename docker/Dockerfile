# By default build with latest from mathworks/matlab (make sure you have docker pull'd this) - but you can over-ride
# this (and the github workflows frequently do) so as to build with different versions of MATLAB
ARG MATLAB_DOCKER_TAG=latest

FROM mathworks/matlab:${MATLAB_DOCKER_TAG}

# By default build from the source repository - if you happen to have forked this source respository then you should 
# ensure you have at least released one version from your repo by tagging and commiting the tag.
ARG SRC_REPO=UniprJRC/FSDA

# This argument MUST be specified in a build of this Dockerfile. See https://github.com/UniprJRC/FSDA/tags for
# valid tags
ARG FSDA_RELEASE

USER root
WORKDIR /root

RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && \
    apt-get install --no-install-recommends --yes wget && \
    apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/matlab
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \
    chmod +x mpm && \
    MATLAB_RELEASE=`ls | grep R*` && \
    ./mpm install --destination=/opt/matlab/${MATLAB_RELEASE}/ --release=${MATLAB_RELEASE} \
        Statistics_and_Machine_Learning_Toolbox \
        Parallel_Computing_Toolbox \
        MATLAB_Coder \
        Optimization_Toolbox  || \
    (cat /tmp/mathworks_root.log && false) && \
    rm -f mpm /tmp/mathworks_root.log

WORKDIR /opt/fsda
# It is critical that the FSDA folder is writable by the user matlab so that the builddocsearchdb function 
# that is run during MATLAB startup can succeed. In additon we would like to know the actual folder name
# for the FSDA toolbox later so write this to the well-known file /opt/fsda/fsda-location.txt
RUN wget -q https://github.com/${SRC_REPO}/archive/refs/tags/${FSDA_RELEASE}.tar.gz && \
    tar -xzf ${FSDA_RELEASE}.tar.gz  && \
    rm ${FSDA_RELEASE}.tar.gz && \
    cd FSDA-${FSDA_RELEASE} && \
    rm -rf .github circleci _* && \
    echo -n `pwd` > ../fsda-location.txt && \
    chown -R matlab .
WORKDIR /

# Seems to be needed for -browser version to install addons - TODO not sure why
RUN mkdir /.Add-Ons && \
    chown -R matlab /.Add-Ons

# One time startup.m that undertakes installation etc. of FSDA in MATLAB
COPY startup.m /home/matlab/Documents/MATLAB/startup.m

USER matlab
WORKDIR /home/matlab