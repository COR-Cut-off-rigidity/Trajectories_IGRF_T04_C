## Precision testing

Test executable disables T04 and IGRF models and sets uniform magnetic field based on settings in `CMakeList.txt` file. Program then performs calculation of the particle's trajectory. Based on the fact that we know precisely how a uniform magnetic field should affect the trajectory of the particle, it is possible to test whether there is a cumulative calculation error in trajectory calculation.

In the current version, there is just a simple check, which observes changes of a particle's trajectory radius. If this change is greater than the maximum allowed (in `CMakeList.txt` file), the test fails.

### Building executable for testing of code precision

1. In `./src/CMakeList.txt` file change `TRAJECTORY_TEST` to `TRUE`;
    - You can specify whether you want to log trajectory points of a particle by setting `PRINT_TRAJECTORY` to `TRUE` or `FALSE`;
    - You can also specify components of a magnetic field by modifying `BX`, `BY` and `BZ` values;
    - `MAX_RD` stands for *maximal radius difference*, which sets the "failure point" of a test.
5. Run `cmake -S ./src -B ./build` command from the root directory of the code;
6. Navigate to build directory: `cd build`;
7. Type `make` to compile the code. Executable binary will be created in `build` directory.


### Outfile (test) example:

```
Particle's rigidity: 0.332
Intensity of magnetic field: X = 0.000000, Y = 86.926332, Z = 0.000000

x y z r theta phi time
12742400.000000 0.000000 -12511.525391 2.000001 1.571778 0.000000 0.000000
12742388.000000 0.000000 -25022.488281 2.000002 1.572760 0.000000 0.000125
12742363.000000 0.000000 -37533.433594 2.000003 1.573742 0.000000 0.000250
```