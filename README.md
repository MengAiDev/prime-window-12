# prime-window-12

High-performance search for a natural number `N` such that the closed interval
`[N, N+2003]` contains **exactly 12 primes**, with

```
N in [1011128158584751, 10098097238186292]
```

## Algorithm

- Deterministic segmented sieve (odds-only bitset, incremental per-block offsets,
  no per-prime division in the hot loop) over the full range.
- Deterministic Miller-Rabin (7 bases, correct for all `u64`) for final verification.
- Event-driven window counter: `f(N) = #{primes in [N, N+2003]}` updated only at
  enter/exit events, so every `N` is checked without per-`N` work.
- 2 worker threads per process, block size 2e7 (empirically fastest).

Throughput ~ 0.8-1.0 G positions/s per 2-core machine.

## Usage

```
g++ -O3 -march=native -funroll-loops -pthread -o scan_range scan_range.cpp
./scan_range <start> <end> [block_size]
```

`scan_range` exhaustively scans `[start, end]` and prints the first `N` found
(or `NO SOLUTION FOUND in range`).  A solution, if printed, is verified by
independent Miller-Rabin before output.

## Parallel computation on GitHub Actions

`.github/workflows/scan.yml` splits the search range into chunks and runs them
in parallel as a matrix of jobs (each on a 2-core runner).  Trigger it via
Actions -> parallel-scan -> Run workflow, with inputs:

- `range_start` / `range_end`: sub-range to cover
- `chunk_size`: positions per job (keep <= ~1.6e13 to finish within the 6 h job limit)
- `max_chunks`: number of parallel jobs (<= 256)

Each job uploads its log as an artifact; a final `collect` job aggregates the
verdicts into `results/summary.txt` on the `results` branch.

`.github/workflows/ci.yml` compiles and smoke-tests the program on every push.

## Math note

`[N, N+2003]` has 12 primes only when 13 consecutive prime gaps (14 primes) sum
to >= 2005.  That is an extreme-tail event at this scale, so solutions (if any)
are very rare; the search is designed for exhaustive, parallel coverage.
