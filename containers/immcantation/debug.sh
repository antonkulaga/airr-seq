#!/bin/bash
docker run  -v /data:/data -v `pwd`/src:/home/magus/notebooks --network=host -i -t quay.io/comp-bio-aging/immcantation jupyter lab .
#docker run --network=host -i -t immcantation/lab:devel
