cmake_minimum_required(VERSION 2.8.10)
project(T04_C C)

# =============================
# BUILD SETTINGS
# =============================

# Trajectory testing (following settings are ignored if TRAJECTORY_TEST is set to FALSE)
set(TRAJECTORY_TEST FALSE)
set(PRINT_TRAJECTORY TRUE)
set(BX 0)
set(BY 86.9263319175807497238)
set(BZ 0)
set(MAX_RD 1)

# TU angle
set(TU 0.01f)

# =============================

set(CMAKE_C_FLAGS "-Wall -Werror -O2") 

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

if(TRAJECTORY_TEST)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-variable -Wno-unused-but-set-variable -DTRAJ_TEST")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_BX=${BX}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_BY=${BY}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_BZ=${BZ}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_MAX_RD=${MAX_RD}")
    if(PRINT_TRAJECTORY)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPRINT_TRAJECTORY")
    endif()
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEF_TU=${TU}")

add_executable(T04_C main.c birk_tot.c birk_tot.h btot.c btot.h deformed.c deformed.h full_rc.c full_rc.h
                            positions.c positions.h recalc.c recalc.h utility.c utility.h test.c test.h)

target_link_libraries(T04_C m)