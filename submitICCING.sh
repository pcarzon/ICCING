#!/bin/bash

module load gcc
module load cmake
module load boost/.1.71.0

cd /projects/jnorhos/pcarzon/ICCING
./iccing $1
