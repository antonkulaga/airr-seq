FROM quay.io/comp-bio-aging/base:latest

WORKDIR /opt
RUN git clone git://github.com/csi-greifflab/Absolut
RUN apt install -y libqt5svg5 libqt5printsupport5 libgsl-dev
RUN apt install -y freeglut3-dev
RUN DEBIAN_FRONTEND=noninteractive apt install -y qtbase5-dev qtcreator
WORKDIR /opt/Absolut/src/bin
RUN qmake ../Absolut/Absolut.pro	#qmake creates a new Makefile embedding Qt libraries locations
RUN make	#This creates Absolut
ENV PATH="/opt/Absolut/src/bin:${PATH}"
WORKDIR /data

CMD ["Absolut"]