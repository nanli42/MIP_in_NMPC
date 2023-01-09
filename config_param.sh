#!/bin/bash

delta_s=0.06
step_N=$1
lookahead=$(echo "$delta_s*$step_N" | bc)

# change step_N and lookahead in source file for generation
cmd_string="%s/OCP ocp( s_start, s_end, .* );/OCP ocp( s_start, s_end, ${step_N} );/g |"
cmd_string+="%s/const double s_end = .*;/const double s_end = ${lookahead};/g |"
cmd_string+="wq";
vim -c "$cmd_string" ../acado_code_generation/src/formulation.cpp

# generate ACADO framework
cd ../acado_code_generation
rm -rf build; mkdir -p build; cd build; cmake ..; make
./main
cd ../../build

# change step_N, step_s etc in source files
Ntot=$(echo "23*$step_N" | bc)
cmd_string="%s/#define Nstep .*/#define Nstep ${step_N}/g |"
cmd_string+="%s/#define Ntot .*/#define Ntot ${Ntot}/g |"
cmd_string+="%s/#define Nobs .*/#define Nobs ${step_N}/g |"
cmd_string+="%s/#define step_s .*/#define step_s ${delta_s}/g |"
cmd_string+="wq";
vim -c "$cmd_string" ../include/step_param.h

vim -c "%s/acado_solve/acado_solve_qpoases/g | wq" ../external/acado_generated/acado_qpoases_interface.hpp
vim -c "%s/acado_solve/acado_solve_qpoases/g | wq" ../external/acado_generated/acado_qpoases_interface.cpp

mkdir -p ../result/track1/N15
mkdir -p ../result/track1/N30
mkdir -p ../result/track2/N15
mkdir -p ../result/track2/N30

bash ../select_mode.sh $2
