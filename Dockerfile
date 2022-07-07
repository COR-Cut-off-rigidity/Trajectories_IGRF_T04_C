FROM debian:11

WORKDIR /usr/src/Trajectories_IGRF_T04_C

COPY . .

RUN apt update && \
    apt install -y gcc cmake make
RUN adduser --no-create-home --disabled-password --gecos "" cor
RUN ./build.sh

RUN mkdir /var/data
VOLUME /var/data
USER cor

ENTRYPOINT [ "./docker-entrypoint.sh" ]