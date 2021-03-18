#/bin/bash
docker run --network=host -i -t quay.io/comp-bio-aging/immcantation jupyter lab .
#docker run --network=host -i -t immcantation/lab:devel