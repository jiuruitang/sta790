FROM rocker/tidyverse:3.5.2

RUN install2.r --deps TRUE\
    mvtnorm \
    RecordLinkage \
    && mkdir /work

WORKDIR /work

ENTRYPOINT ["Rscript"]