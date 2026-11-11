# Assembly-guided genotype packing — 2026-11-10

This measures an incremental change from **68dade2**, not another comparison with 1.x. The numerical formulas, logarithms, LD solver, pair-counting dispatch thresholds, NEON counters and AVX2/AVX-512 counters are unchanged.

## Change and assembly evidence

The old packing loop selected one of three category planes, branched on every genotype and performed a read/modify/write of a bit-plane word for every observed sample. Packing now constructs all three words in registers and writes each word once. Full blocks use 16-sample NEON or baseline x86-64 SSE2 comparisons; partial blocks and other architectures use exact integer comparisons and shifts. Reusing allocated packed storage no longer clears words that will all be overwritten.

On ARM64 the packing probe emits `cmeq`, `and`/`orr`, `addv` and scalar shifts to assemble exact sample-position masks. The full-word path has register bit masks and sequential input loads, with no per-sample scatter or data-dependent category branch. Function prologue register saves and scalar-tail handling still exist; this is not a claim that the entire function has no stack accesses. On x86-64, `pcmpeqd`/`movmskps` gather the masks using SSE2, which is available in the architecture baseline. No AVX requirement is added to packing. Counting continues to select POPCNT, AVX2 or AVX-512 independently as before.

The reproducible [assembly probes](../../assembly/README.md) now expose packing too. Compiler details, emitted instruction counts and executable SHA-256 values are in [metadata.json](metadata.json).

## Measured ARM64 results

macOS arm64, Homebrew Clang 22.1.8, `-O3`, C++17, seed 2026. Engine measurements use the existing Release/OpenMP build with **one thread**. Synthetic observed genotypes are uniformly distributed over 0/1/2. Input generation and I/O are excluded; engine timing includes summaries, packing, caches, pair calculations and matrix allocation. These are not file-to-output CLI timings.

| Workload | Packing speedup | MI engine | RMI engine | r² engine |
| --- | ---: | ---: | ---: | ---: |
| 4,096 samples × 128 variants, 25% missing | 6.30× | 2.03× | 1.70× | 1.69× |
| 65,536 samples × 128 variants, 25% missing | 6.94× | 2.27× | 2.19× | 2.25× |
| 512 samples × 768 variants, 25% missing | See separate packing dimensions | 1.09× | 1.04× | 1.02× |

Packing and engine dimensions differ in the small-sample tests: the packing-only test uses 384 variants, while the engine uses 768. The tables must not be combined as if they were the same workload. Full measured results, including zero/100% missing packing cases, are in [packing-arm64.tsv](packing-arm64.tsv) and [engine-arm64.tsv](engine-arm64.tsv).

Packing measurements take one warmup and 15 timed builds per invocation, then the median of five alternating before/after invocations. Engine measurements use the existing executable's one warmup and three timed evaluations for each method, then the median of three alternating before/after invocations. Repeated packing builds reuse vector capacity; engine evaluations construct a new engine each time. Timings below a millisecond are sensitive to host scheduling and CPU frequency. The larger workloads provide stronger evidence than small percentage changes in short runs.

Every MI, RMI and r² output matrix in the six engine workloads was compared byte for byte, including NaN entries, and matched the baseline. Independent LD checks and regression tests cover D/D′ separately; this table does not claim new D/D′ timings.

## x86 verification and limits

The updated bit-packing and counting tests also pass in macOS x86-64 execution through Rosetta, including every invalid-code lane and every output bit. [packing-rosetta.tsv](packing-rosetta.tsv) records supplementary packing timings. These are translation timings, **not native x86 hardware measurements**. Ordinary mixed data is faster in that run; the all-missing case is slower because the old scalar branch skips essentially all work. There is no universal speedup claim for every input distribution.

Native ARM64 GCC 16 and Clang builds pass the tests. AddressSanitizer/UndefinedBehaviorSanitizer, all six OpenMP CTest groups, all five serial CTest groups, all six R regression scripts with CLI/VCF integration and R CMD check pass. The R package and local CLI were reinstalled. AVX2/AVX-512 assembly checks still pass; this session did not execute those vector counters on supporting physical x86 hardware.

## Reproduce

Use the same compiler and flags for both snapshots. From the repository root:

```bash
mkdir -p build/packing-before build/packing-after
# Extract headers only; this does not change the worktree.
git archive 68dade2 GWLD-R/inst/include/gwld | tar -x -C build/packing-before
CXX=/opt/homebrew/opt/llvm/bin/clang++
"$CXX" -std=c++17 -D_REENTRANT -O3 \
  -I build/packing-before/GWLD-R/inst/include \
  tests/assembly/packing_benchmark.cpp -o build/packing-before/packing_benchmark
"$CXX" -std=c++17 -D_REENTRANT -O3 -I GWLD-R/inst/include \
  tests/assembly/packing_benchmark.cpp -o build/packing-after/packing_benchmark
```

Compile the existing engine benchmark against each header snapshot. These commands match the LLVM/OpenMP setup used on this macOS host:

```bash
LLVM_PREFIX=/opt/homebrew/opt/llvm
ARMADILLO_INCLUDE=/opt/homebrew/include
for revision in before after; do
  if [ "$revision" = before ]; then
    SHARED_INCLUDE=build/packing-before/GWLD-R/inst/include
  else
    SHARED_INCLUDE=GWLD-R/inst/include
  fi
  "$CXX" -std=c++17 -O3 -DNDEBUG -D_REENTRANT \
    -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS \
    -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU \
    -fopenmp=libomp -I "$ARMADILLO_INCLUDE" -I "$SHARED_INCLUDE" \
    GWLD-C++/tests/kernel_benchmark.cpp \
    -L "$LLVM_PREFIX/lib" -Wl,-rpath,"$LLVM_PREFIX/lib" \
    -o "build/packing-$revision/gwld_kernel_benchmark"
done
```

Keep the compiler, Armadillo and OpenMP libraries identical across builds. For other platforms, use their native dependency locations and OpenMP flags. The production CMake build's `flags.make` and `link.txt` also record its platform-specific options.

```bash
python3 tests/assembly/compare_packing.py \
  --before build/packing-before --after build/packing-after \
  --output build/packing-comparison
```

The driver repeats the measurements above and checks byte-identical engine matrices. Run it without competing builds/tests. It writes median timings and the last pair of compared matrices. For other operating systems, use native compiler/library locations and the same optimization settings in both builds; do not compare Rosetta time directly with native ARM64 time.

## Follow-up review

The implementation, table arithmetic and recorded executable hashes were rechecked on 2026-11-10. OpenMP CTest passed 6/6 groups and a freshly rebuilt serial configuration passed 5/5. The checked-in comparison driver was executed against the recorded before/after binaries; all six workloads again produced byte-identical MI, RMI and r² matrices. The documented standalone engine compile command was smoke-tested. The historical timings above are retained rather than replaced with every repeated run.

The full R example, paper-subset workflow and serial R installation were validated in the earlier review, not rerun in this packing follow-up. Current packing-specific R evidence consists of the six OpenMP-package scripts and the clean R CMD check described above.
