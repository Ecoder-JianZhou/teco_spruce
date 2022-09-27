FROM centos:7
MAINTAINER Mark Stacy <markstacy@ou.edu>

RUN mkdir /source
COPY . /source
WORKDIR /source
RUN yum install -y gcc gcc-gfortran
RUN gfortran -o /source/teco_spruce /source/TECO_SPRUCE.f90
RUN chmod +x /source/teco_spruce
ENTRYPOINT ["./teco_spruce"]
CMD ["./input/SPRUCE_pars.txt", "./input/SPRUCE_forcing.txt","./input/SPRUCE_obs.txt","/source/output/","0"]



