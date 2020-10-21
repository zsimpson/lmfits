FROM ubuntu:latest

RUN apt-get update && apt-get install -y build-essential libblas-dev liblapack-dev f2c libgfortran5

COPY ./levmar-2.6 /app/levmar-2.6
RUN cd /app/levmar-2.6 && make

COPY . /app
RUN cd /app && gcc gauss2.c -I./levmar-2.6 -L./levmar-2.6 -llevmar -lm -llapack -lblas -o lmfitter
# -O3 for optimzation

WORKDIR /app
ENTRYPOINT /bin/bash
