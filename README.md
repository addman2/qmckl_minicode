# qmckl_minicode

Install with openmp offload support

```bash
#gfortran
cmake -S <source> -B <build> -DOFFLOAD_FLAGS="-fopenmp -foffload-options='-lm -lgfortran ...'"

#nvfortan
cmake -S <source> -B <build> -DOFFLOAD_FLAGS="-mp=gpu -gpu=cuda11.4"

# ... you know the drill
```
