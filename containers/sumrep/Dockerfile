FROM quay.io/comp-bio-aging/rstudio

MAINTAINER antonkulaga@gmail.com

WORKDIR /opt
RUN git clone https://github.com/matsengrp/sumrep

ADD install.R /tmp/
RUN R -f /tmp/install.R

#install extra packages
RUN pip install ./pygor

WORKDIR /home/rstudio