#! /bin/bash
docker build -f Dockerfile . -t methylfastq && \
docker run -it --rm --entrypoint /bin/bash -w /tool -v "$(pwd):/tool" -v ~/.aws:/root/.aws methylfastq
