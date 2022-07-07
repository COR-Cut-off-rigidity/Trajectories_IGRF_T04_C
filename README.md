# Trajectory simulation of cosmic rays in the Earth's magnetosphere

This program produces calculated values of cut-off rigidities for cosmic ray trajectories, which describes spectrum of allowed and forbidden rigidities.

Program uses **Tsyganenko 04** model for simulating Earth's external geomagnetic field and **IGRF** (gen. 9 - 13) for simulating Earth's internal geomagnetic field. It was developed as refactoring of *Institute of Experimental Physics, Slovak Academy of Sciences* Fortran code.

> This is part of [cor.crmodels.org](https://cor.crmodels.org/) project.

## Build instructions

For building your own executable binary, you need to have installed the following prerequisites:

- [CMake](https://cmake.org/) - v2.8.10 or newer
- [Make](https://www.gnu.org/software/make/) - v4.2.1 or newer
- [GCC](https://gcc.gnu.org/) - v9.3 (recommended)

Assuming you have already cloned this repository (or downloaded and extracted ZIP archive):

1. Navigate to the project directory
2. Type `./build.sh` into the teminal
3. Executable binary `Trajectories_IGRF_T04_C` will be created inside `build` directory

Alternatively you can build a docker image and run this simulation in a containerized environment. In this case, you only need to have [Docker](https://www.docker.com/get-started/) installed on your machine.

1. Navigate to the project directory
2. Type `./build-image.sh` into the teminal. This script executes `docker build` command in order to build the image called `trajectories_igrf_t04_c` from the source with the `latest` tag.

## How to run the simulation

This program requires several input parameters:

1. `<infile>` - file, that contains data necessary for the simulation,
2. `<outfile>` - file, into which the results will be saved,
3. `<igrf_ver>` - Version of IGRF coefficients (value: 9 - 13),
4. `<seq/par>` - Sequential (seq) or parallel (par) computation,
5. `<number of steps>` - Max. number of steps for single trajectory calculation.

Assuming you already have an executable binary, you can execute the program by running:
`./Trajectories_IGRF_T04_C <infile> <outfile> <igrf_ver> <seq/par> <number of steps>`.

More specific command for running the simulation with the input file located in the same directory as the executable: `./Trajectories_IGRF_T04_C infile outfile 13 par 25000`. In this case, 13th generation of the IGRF model is used, simulation will run in parallel mode with maximum of 25000 steps for single trajectory calculation.

## How to run the simulation in containerized environment

If you want to run the simulation in docker container, it is necessary to create bind mount or named volume for data persistance. We recommend using bind mount as it is easier to work with. For that, you need to create directory called `data`, which will be mapped to container's file system. In this directory you need to have appropriate input file created.

Assuming you have already built the image, you can use following command to create the container and run the simulation: 

```sh
docker run --name trajectories_igrf_t04_c -v $(pwd)/data:/var/data -e INFILE=infile trajectories_igrf_t04_c
```

This command will create container named `trajectories_igrf_t04_c`, from the image `trajectories_igrf_t04_c`, with the directory `./data` mapped to the `/var/data` directory inside container and run the simulation. When you have the container created, you can run the simulation again by executing: `docker start -a trajectories_igrf_t04_c`.

In order to remove existing container, run `docker rm trajectories_igrf_t04_c`. There is a flag `--rm` you can use with `docker run` command. This will cause the created container to 
be automatically removed after it exits.

Environment variables:

- `INFILE` **(required)** - Name of the infile
- `OUTFILE` - Name of the outfile (default: `outfile`)
- `IGRF` - Version of IGRF coefficients (default: `13`)
- `MODE` - Mode of the computation. Use `seq` for sequential and `par` for parallel calculation (default: `par`)
- `STEPS` - Max. number of steps for single trajectory calculation (default: `25000`)

## Input file example

```
  2.5000  -1.   20.0000
  1.00     48.66     20.53
            48.66     20.53
 2000  3 28  88 16 00 00
   100   1   1   0.10
   -2   0.37   2.30  -0.80
   0.04   0.03   0.06   0.01   0.03   0.02
-1.00
```

1. Line: starting rigidity; type of particle (-1 for proton); ending rigidity
2. Line: radius in `Re` - (Earth's radius); geographic latitude and longitude (of the starting point of the trajectory)
3. Line: geographic latitude and longitude for the incoming direction of the particle
4. Line: year; month; day; day in year; hours; minutes; seconds
5. Line: First 3 numbers - parameters for trajectory division; rigidity step
6. Line: `Dst` index (in `nT`); `pdyn` - dynamic pressure of solar wind (in `nPa`) at given date and time; intensity of `y` and `z` component of interplanetary field at given date and time
7. Line: W1 - W6 input parameters for external geomagnetic field model (Tsyganenko) which describes prehistory of geomagnetic field.
8. Line: `-1.00` marks the end of the input file.

## Output file example

```
ASYMPTOTIC COORDINATES
calculated by model of exter.field T05
 Station with geo.latitude: 48.660 & longitude: 20.530 & radius: 1.00000
 Direction of trajectory with latitude: 48.660 & longitude: 20.530
 Date: 2000 3 28 time: 16 hod 0 min 0 sec
 Starting rigidity: 2.5 GV Epsilon=0.100000
 Limit of total number of steps: 25000

 rig : v : rad : eth : efi : ath : afi : time : length
  3.500000 0.9659288526 24.348516 1.976 177.357 -1.295 193.369 1.397261 406080.81
  3.900000 0.9722833633 11.134421 -14.778 329.334 -10.356 3.128 0.678627 197938.64
  ...
  19.900000 0.9988912344 23.895048 14.815 66.211 12.526 69.865 0.483876 148820.27
  20.000000 0.9989024401 23.824049 15.035 66.083 12.740 69.722 0.475542 147421.17
  CUTOFF with rigidities P(S),P(C),P(M) are:
3.50000 3.90000 3.80000
```

This file consists of 2 parts:
- Header - contains simulation metadata (geographic coordinates, date and time, Epsilon - step in rigidity, maximum number of steps used for backtracing of one trajectory)
- Calculated data - contains trajectory parameters of cosmic rays particle with given rigidity in given time and location. Last line of the outfile includes 3 values, which describes lower, upper and effective cutoff rigidity.

## Program flowchart

![Program flowchart](./docs/flowchart.jpg)
