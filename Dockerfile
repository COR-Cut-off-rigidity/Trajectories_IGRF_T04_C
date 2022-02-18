FROM debian:latest

WORKDIR /usr/src/Trajectories_IGRF_T04_C

COPY . .

RUN apt update && apt install -y gcc cmake make
RUN ./build.sh