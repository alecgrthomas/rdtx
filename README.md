# rdtx
RDTX code - calculates radiation from charged particles. 

## Installation instructions

Makefile is written for GCC, tested with 10.0.1 on MacOS

To build
`make`

To run
`./bin/rdtx_beta <deck>`
where `<deck>` is the name of the input deck containing the parameters you want to run.

If `<deck>` is blank / unspecified it will look for `defaults.dat`
