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

### 2. Coordinate text silently rounded through floating point

The old parser used `std::stold`. On the reviewed ARM64 macOS toolchain, `long double` has the same integer precision as `double`. Reproductions before the fix were:

| Input text | Old result | Correct handling |
| --- | --- | --- |
| `9007199254740993` | `9007199254740992` | Preserve the exact integer. |
| `1.00000000000000001` | `1` | Reject a fractional coordinate. |
| `4.611686018427387903e18` | Rejected after rounding | Accept the exact safe maximum `4611686018427387903`. |

Ordinary integers now use integer `from_chars`. Decimal/scientific notation is parsed as an exact decimal coefficient and exponent; fractional results and range overflow are rejected. The existing nonnegative coordinate range, capped at `INT64_MAX/2`, remains in effect. Numeric metadata such as `5e9` remains accepted. Hexadecimal floating-point spellings are rejected; coordinates use decimal notation.

The core tests include 1,000 deterministic integer/scientific round trips and rejection of fractional suffixes. CLI tests check the actual PLINK-to-genotype export. R tests pass character metadata above 2^53 through pair calculation and verify an exact one-base distance. This does not recover precision already lost when a caller constructs a large coordinate as an R double; exact high-range metadata must be supplied as text. VCF, PED and BED readers now switch the whole POS column to exact decimal text when any coordinate exceeds 2^53; ordinary coordinates retain numeric storage. Invalid fractional MAP/BIM positions are rejected at read time.

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
- R package installation and all six R test scripts in OpenMP and serial builds, including CLI/R VCF integration and exact-coordinate reader checks.
- R CMD check with no errors, warnings or notes.
- The complete R example, paper-subset validation and native interruption tests.

This review does not execute AVX2 or AVX-512 on native supporting hardware. Docker's daemon was unavailable for a fresh QEMU run. Earlier platform measurements are documented separately and do not replace runtime validation of the revised feature checks on every x86 configuration. The simulated feature-selection tests and generated x86 assembly cover the specific dispatch defect found here.
