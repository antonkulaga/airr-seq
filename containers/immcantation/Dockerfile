FROM immcantation/suite:devel

ENV CRAN_URL https://cloud.r-project.org/
#ENV DISABLE_AUTH=true
ENV PASSWORD=rstudio

ADD https://s3.amazonaws.com/rstudio-server/current.ver /tmp/ver

RUN yum install rstudio-server -y
USER rstudio
WORKDIR /home/rstudio

EXPOSE 8787

#ENTRYPOINT ["/usr/bin/rserver"]
CMD "/usr/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0"
