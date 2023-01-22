
Results obtained by

```shell
hyperfine "cargo run --bin rejection_sampling_sir_mf" --warmup 3
```

Benchmark 1: cargo run --bin rejection_sampling_sir_mf
  no buffered writer
  Time (mean ± σ):     19.967 s ±  2.075 s    [User: 3.942 s, System: 15.880 s]
  Range (min … max):   18.398 s … 25.249 s    10 runs

Benchmark 1: cargo run --bin rejection_sampling_sir_mf
  default buffered writer
  Time (mean ± σ):      1.226 s ±  0.333 s    [User: 0.890 s, System: 0.361 s]
  Range (min … max):    0.776 s …  1.765 s    10 runs

Benchmark 1: cargo run --bin rejection_sampling_sir_mf
  buffered with 100MB
  Time (mean ± σ):     864.2 ms ± 232.7 ms    [User: 643.4 ms, System: 222.5 ms]
  Range (min … max):   731.7 ms … 1347.4 ms    10 runs
