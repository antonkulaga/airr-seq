#/bin/bash
docker run --network=host -v /data:/data -i -t quay.io/comp-bio-aging/immcantation jupyter lab .