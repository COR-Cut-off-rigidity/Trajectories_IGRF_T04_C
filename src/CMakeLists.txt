cmake_minimum_required(VERSION 2.8.10)
project(Trajectories_IGRF_T04_C C)

# =============================
# MODULES
# =============================

include(CMakeDependentOption)

# =============================
# BUILD OPTIONS
# =============================

option(TRAJECTORY_TEST
    "Build test executable" OFF
)

cmake_dependent_option(CHECK_R
    "Enable radius value check during runtime" OFF
"TRAJECTORY_TEST" OFF)

cmake_dependent_option(WRITE_TRAJECTORY
    "Write trajectory values to outfile" ON
"TRAJECTORY_TEST" OFF)

# =============================
# BUILD VARIABLES
# =============================

# These values will be written to outfile (if PRINT_TRAJECTORY
# option is enabled). Possible values: X;Y;Z;R;TH;F;TIME;MAX_RD
set(WRITE_VALUES "X;Y;Z;MAX_RD")

# This value specifies number of trajectory steps between
# writes to outfile (if PRINT_TRAJECTORY option is enabled)
set(WRITE_STEP 100)

# Values of magnetic field vector components used in test build
set(BX 0)
set(BY 86.9263319175807497238)
set(BZ 0)

# This value specifies max. difference of current radius and radius
# specified in infile (in meters). This value is used when CHECK_R option is enabled. 655000
set(MAX_RD 655000)

# This value specifies value of angle "epsilon". This value is ignored
# when TRAJECTORY_TEST option is enabled, but it is possible to specify
# this value as a command line argument
set(TU 0.01f)

# =============================

set(CMAKE_C_FLAGS "-Wall -Werror -O2") # -mavx2 -mavx -march=native

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

if(TRAJECTORY_TEST)
    message(STATUS "TEST MODE IS ENABLED!!!")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-variable -Wno-unused-but-set-variable -DTRAJ_TEST")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_BX=${BX} -DDEF_BY=${BY} -DDEF_BZ=${BZ}")

    if(CHECK_R)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DCHECK_R -DDEF_MAX_RD=${MAX_RD}")
    endif()

    if(WRITE_TRAJECTORY)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWRITE_TRAJ -DDEF_WRITE_STEP=${WRITE_STEP}")
        foreach(VAL IN LISTS WRITE_VALUES)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWRITE_VAL_${VAL}")
        endforeach()
    endif()
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_TU=${TU}")

add_executable(Trajectories_IGRF_T04_C main.c birk_tot.c birk_tot.h btot.c btot.h deformed.c deformed.h full_rc.c full_rc.h
                            positions.c positions.h recalc.c recalc.h utility.c utility.h test.c test.h)

target_link_libraries(Trajectories_IGRF_T04_C m)

# =============================
# CLEANUP
# =============================

unset(TRAJECTORY_TEST CACHE)
unset(CHECK_R CACHE)
unset(WRITE_TRAJECTORY CACHE)