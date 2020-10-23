FROM ubuntu:latest

RUN apt-get update && apt-get install -y build-essential libblas-dev liblapack-dev f2c libgfortran5

COPY ./levmar-2.6 /app/levmar-2.6
RUN cd /app/levmar-2.6 && ENV_OPTS="-fPIC" make

COPY . /app
RUN cd /app && LEVMAR_FOLDER=/app/levmar-2.6 DST_FOLDER=/app ./build.sh

WORKDIR /app
ENTRYPOINT /bin/bash
