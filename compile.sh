#!/usr/bin/env bash
#
export BUILD_DIR=build-docker
export CUDA_VISIBLE_DEVICES=0

rm -rf ${BUILD_DIR}
cmake -S . -B ${BUILD_DIR} -DCMAKE_C_COMPILER=nvc \
                           -DCMAKE_Fortran_COMPILER=nvfortran \
                           -DOFFLOAD_FLAGS="-mp=gpu -gpu=cc75" \
                           -DOLD_QMCKL_GPU_INTERFACE=ON \

cmake --build ${BUILD_DIR} --verbose
#ctest --test-dir build -VV -R "ComparatorWithoutPseudo"
#ctest --test-dir ${BUILD_DIR} -VV
${BUILD_DIR}/comparator
