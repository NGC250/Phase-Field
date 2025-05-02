#!/bin/zsh

RUN_DIR="$(pwd)"
SOURCE_DIR="$RUN_DIR/../src"
PYTHON_DIR="$RUN_DIR/../Python"

OMP_DIR="$(brew --prefix libomp)"
FFTW_DIR="$(brew --prefix fftw)"

clang -Wall -Wextra -std=c11 -g -o domain "$SOURCE_DIR/main.c" -L"$OMP_DIR/lib" -L"$FFTW_DIR/lib" -I"$OMP_DIR/include" -I"$FFTW_DIR/include" -lfftw3 -lm -fopenmp #-lfftw3_threads

chmod +x domain

./domain

source "$PYTHON_DIR/analysis/bin/activate"

python "$PYTHON_DIR/scripts/analysis.py"

deactivate
