# Assembly-assisted code review — 2026-11-10

This review starts from `bb2a303` and covers production numerical headers, packed counting and dispatch, native R adapters, genotype encoding, input/output and the R workflows. It combines source inspection, optimized assembly, executable regression tests and independent LD validation. It is not a formal proof that every possible input or compiler configuration is correct.

## Reproduce the assembly inspection

`audit_kernels.py` compiles small, named probes against the actual shared headers. It requires Python 3 and Clang or GCC. It saves the compiler command, version, assembly and instruction counts; it does not alter the production compiler flags.

```bash
# Native architecture; choose the compiler explicitly for reproducibility.
python3 tests/assembly/audit_kernels.py \
  --cxx /opt/homebrew/opt/llvm/bin/clang++ --output build/assembly/native

# macOS cross-assembly only: does not execute AVX2 or AVX-512 instructions.
python3 tests/assembly/audit_kernels.py \
  --cxx /opt/homebrew/opt/llvm/bin/clang++ --arch x86_64 \
  --output build/assembly/x86
```

On macOS the script obtains the SDK with `xcrun`. On Linux use a native compiler through `--cxx`; the macOS `--arch` option is not a Linux cross-toolchain configuration. The probes are inspection-only and are not linked into the CLI or R package. MSVC is not supported by this compiler-driver script.

Generated files are `kernels.s`, `instructions.json` and `compiler.log`. Recorded instruction counts are in [recorded/2026-11-10.json](recorded/2026-11-10.json). Checks confirm the expected NEON, AVX2 byte-shuffle and AVX-512 vector-popcount instructions when those functions are emitted. On x86 the portable fallback is also checked for unintended POPCNT instructions. A passing instruction check establishes emitted code structure, not runtime correctness on an unavailable CPU.

## Findings and fixes

### 1. Scalar POPCNT hidden inside vector kernels

Apple Clang 17 emitted scalar `popcntq` instructions in the tails of functions annotated only with `target("avx2")` or `target("avx512f,avx512vpopcntdq")`. The source tail used portable SWAR arithmetic, but the optimizer recognized that idiom. Full and complete-pair variants were affected. The old dispatcher tested the vector features without testing scalar POPCNT. LLVM also emitted AVX2 reduction/tail instructions inside AVX-512 functions.

Ordinary AVX2-capable machines normally provide scalar POPCNT too, but the dispatcher must match the actual instructions: a feature can be masked independently in a virtual machine. The revised functions explicitly declare `popcnt`, and dispatch requires it before either vector implementation can be selected; AVX-512 selection additionally requires AVX2. Without scalar POPCNT, the portable fallback is used. This avoids depending on compiler-specific feature implications.

`simd_test.cpp` checks the selection result for simulated feature sets where AVX2/AVX-512 are advertised but scalar POPCNT is absent. It compares pointers without calling instructions unsupported by the executing CPU. Actual forced-kernel tests also check the complete required feature set before execution.

### 2. Keep coordinate handling scoped to genomic data

C++ coordinates remain `int64_t`; R readers consistently return numeric POS columns. There is no automatic switch to character columns and no special handling for coordinates above 2^53, which are outside the intended genomic workload.

Ordinary integer text is parsed directly into `int64_t` with `from_chars`. A short fallback accepts scientific notation produced by R for ordinary numeric coordinates. Nonnegative, integral and range checks remain in place. The existing `INT64_MAX/2` bound protects subsequent distance arithmetic.

The initially added arbitrary decimal parser and oversized-coordinate text adapters were removed after scope review. Tests retain ordinary genomic coordinates, scientific notation, invalid fractional/negative input and the existing 5-billion-base example. R reader metadata keeps a stable numeric type.

## Assembly observations

The inspected native build used LLVM 22.1.8, C++17, Release optimization on macOS arm64. x86 assembly was also generated with LLVM 22.1.8 and Apple Clang 17.

| Path | Observation |
| --- | --- |
| Genotype summaries | Integer comparisons and counter reductions are auto-vectorized; LLVM reports width four with interleave count two. |
| Complete-pair NEON loop | `and`, `cnt` and `uadalp`; four counters stay in vector registers. No stack accesses occur in the inspected inner loop. The four-double result uses the ARM64 homogeneous floating-point return ABI. |
| Full NEON counting | Vector population counts and bounded accumulator reductions are present. Tail handling remains separate. |
| AVX2 counting | `vpshufb` implements the nibble lookup; scalar POPCNT appears in tails and is now an explicit dispatch requirement. |
| AVX-512 counting | `vpopcntq` is present; AVX2 and scalar POPCNT are also part of the declared/checkable requirements. |
| x86 call boundary | `vzeroupper` is present in inspected vector routines. Public kernel boundaries pass scalar/pointer arguments and ordinary structs, not ISA-specific vector arguments. |
| RMI/parallel linkage | The actual native CLI imports `lgamma_r` and OpenMP runtime entry points; it does not import `signgam`. |

A missed auto-vectorization remark on a loop containing explicit NEON intrinsics does not mean that the loop is scalar. Conversely, scalar contingency counting retains data-dependent indexed updates; it is used for small inputs where packing overhead is not necessarily worthwhile. No approximate logarithms, fast-math flags or statistical formula changes were introduced.

## Numerical and performance impact

The native ARM64 calculation benchmark executable before and after the fixes had exactly the same SHA-256:

```text
3b344b9c8a31cf1f6c209afb820a7c90f37414395b8c8c356dba99407d6c2b34
```

The coordinate parser and R metadata adapters are unused by that synthetic benchmark, and the x86 dispatch changes are excluded from its ARM64 build. Therefore this review makes no new calculation-speed claim. Existing performance measurements remain separate from these correctness fixes.

The [independent LD oracle](../ld/README.md) still passes all 37 recorded checks. The LD objective and solver are unchanged by this review.

## Verification and limits

Completed checks on 2026-11-10:

- Six CTest groups in the OpenMP build, including independent R LD verification; five groups in the serial build.
- Core and packed-count tests under AddressSanitizer and UndefinedBehaviorSanitizer.
- Native GCC 16 core/SIMD tests.
- macOS x86_64 core and SIMD/dispatch tests executed through Rosetta; only portable/POPCNT runtime paths are available there.
- Native NEON execution and both native/x86 assembly inspections.
- R package installation and all six R test scripts in OpenMP and serial builds, including CLI/R VCF integration and numeric-coordinate reader checks.
- R CMD check with no errors, warnings or notes.
- The complete R example, paper-subset validation and native interruption tests.

This review does not execute AVX2 or AVX-512 on native supporting hardware. Docker's daemon was unavailable for a fresh QEMU run. Earlier platform measurements are documented separately and do not replace runtime validation of the revised feature checks on every x86 configuration. The simulated feature-selection tests and generated x86 assembly cover the specific dispatch defect found here.

## Follow-up packing optimization

The subsequent packing change starts from `68dade2`. It replaces per-sample category scatters with register-built bit-plane words and 16-sample NEON/SSE2 blocks. The probes now include `gwld_pack_probe`; emitted packing compares/movemasks are recorded alongside counting instructions. See [the incremental benchmark report](../benchmarks/2026-11-10/README.md) for timing, byte-level matrix comparisons and architecture limits. The unchanged-binary statement above applies to the earlier dispatch/coordinate review, not this later packing optimization.
