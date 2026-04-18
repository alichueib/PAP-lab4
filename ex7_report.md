# Exercise 7 - Scalability Measurements

## Goal

Exercise 7 asks for scalability measurements of the implemented communication schemes. I prepared a benchmark script that runs exercises 1 to 6 with output disabled, records the `Total time` printed by the program, and writes the results as CSV.

## Method

The script `bench_ex7.sh` supports:

- strong scaling: fixed global problem size, `MODE=strong`;
- weak scaling: problem size grows with process count, `MODE=weak`;
- configurable process counts with `PROCS`;
- configurable exercises with `EXERCISES`;
- repeated runs with `REPEATS`;
- configurable config file with `CONFIG`.

Default command:

```sh
./bench_ex7.sh
```

Example for Grid'5000 or a larger node:

```sh
MODE=strong PROCS="1 2 4 8 16" REPEATS=3 ./bench_ex7.sh
MODE=weak PROCS="1 2 4 8 16" REPEATS=3 OUT=ex7_weak_results.csv ./bench_ex7.sh
```

The program is run with `-n/--no-out`, so timings measure computation and communication without output-file writing.

## Plotting

If `gnuplot` is available, generate the strong-scaling plots with:

```sh
gnuplot plot_ex7.gnuplot
```

This produces:

- `ex7_strong_scaling.pdf`;
- `ex7_speedup.pdf`.

## Local Strong-Scaling Sample

I ran one local strong-scaling sweep with the default `config.txt`, output disabled, and process counts 1, 2, 4, and 8. The results are stored in `ex7_results.csv`.

| Exercise | 1 process | 2 processes | 4 processes | 8 processes | Speedup at 8 processes |
| --- | ---: | ---: | ---: | ---: | ---: |
| ex1 | 6.071420 | 5.101910 | 0.818613 | 0.625944 | 9.70 |
| ex2 | 5.098520 | 3.829790 | 1.438200 | 1.714730 | 2.97 |
| ex3 | 3.685940 | 2.588990 | 0.922357 | 0.623646 | 5.91 |
| ex4 | 2.792380 | 2.644470 | 1.034040 | 1.923160 | 1.45 |
| ex5 | 4.585710 | 2.496780 | 0.885435 | 0.752494 | 6.09 |
| ex6 | 4.609580 | 2.591670 | 1.337920 | 1.838230 | 2.51 |

These numbers are only a local sample with one repetition, so they should not be over-interpreted. A final report should use several repetitions on a dedicated machine or Grid'5000 node and report averages or medians.

In this run, exercise 1, exercise 3, and exercise 5 scale best up to 8 processes. Exercises 2, 4, and 6 improve up to 4 processes but regress at 8 processes, which suggests that communication overhead and local-machine scheduling noise become significant for this problem size.

## Notes for the report

For the final assignment report, include:

- the machine/node type and CPU count;
- the MPI implementation;
- the process counts tested;
- whether the experiment is strong or weak scaling;
- execution-time and speedup graphs;
- a short interpretation comparing exercises 1 to 6.

Expected qualitative behavior:

- exercises 1 to 3 use 1D splitting, so the ghost area grows less favorably as process count increases;
- exercises 4 to 6 use 2D splitting, reducing the ghost-cell surface for larger process counts;
- exercise 5 avoids manual row copy using MPI datatypes;
- exercise 6 uses non-blocking communication, but the benefit depends on MPI progress and whether communication overlaps useful work.
