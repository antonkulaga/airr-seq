#!/bin/bash
docker run  -v /data:/data -v `pwd`/src:/home/magus/notebooks/src --network=host -i -t immcantation_local jupyter lab --ip=0.0.0.0 .
#docker run --network=host -i -t immcantation/lab:devel
