Flags -

for C : -lm
for CUDA : -lcusolver

command C : gcc filename.c -lm
command CUDA : /use/local/cuda-10.0/bin/nvcc filename.cu -lcusolver
