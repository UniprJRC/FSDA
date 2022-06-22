ARG FSDA_RELEASE=8.5.28

FROM mathworks/matlab:r2022a
ARG FSDA_RELEASE
ARG MATLAB_RELEASE=R2022a

USER root
WORKDIR /root

RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && \
    apt-get install --no-install-recommends --yes wget && \
    apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*

RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \
    chmod +x mpm && \
    ./mpm install --destination=/opt/matlab/${MATLAB_RELEASE}/ --release=${MATLAB_RELEASE} \
        Statistics_and_Machine_Learning_Toolbox \
        Parallel_Computing_Toolbox \
        MATLAB_Coder \
        Optimization_Toolbox && \
    rm -f mpm

WORKDIR /opt/fsda
RUN wget -q https://github.com/josmartin/FSDA/archive/refs/tags/${FSDA_RELEASE}.tar.gz && \
    tar -xzf ${FSDA_RELEASE}.tar.gz  && \
    rm ${FSDA_RELEASE}.tar.gz && \
    cd FSDA-${FSDA_RELEASE} && \
    rm -rf .github circleci _* && \
    chown -R matlab .
WORKDIR /


# One time startup.m that undertakes installation etc. of FSDA in MATLAB
COPY startup.m /opt/matlab/${MATLAB_RELEASE}/toolbox/local

# TO REMOVE
COPY license.dat /opt/matlab/${MATLAB_RELEASE}/licenses/license.dat

USER matlab
WORKDIR /home/matlab