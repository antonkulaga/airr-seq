FROM immcantation/suite:devel

ENV CRAN_URL https://cloud.r-project.org/

ADD https://s3.amazonaws.com/rstudio-server/current.ver /tmp/ver

RUN wget https://download2.rstudio.org/server/centos8/x86_64/rstudio-server-rhel-1.4.1103-x86_64.rpm && yum install rstudio-server-rhel-1.4.1103-x86_64.rpm

RUN set -eo pipefail \
      && clir update \
      && useradd -m -d /home/rstudio -g rstudio-server rstudio \
      && echo rstudio:rstudio | chpasswd \
      && echo "r-cran-repos=${CRAN_URL}" >> /etc/rstudio/rsession.conf

EXPOSE 8787

ENTRYPOINT ["/usr/lib/rstudio-server/bin/rserver"]
CMD ["--server-daemonize=0", "--server-app-armor-enabled=0"]