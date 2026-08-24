# Implementation Plan: Ancilla-Assisted Four-T Toffoli Lowering in NWQEC

## Document purpose

This document is an implementation-ready specification for adding an exact, ancilla-assisted, four-\(T\) Toffoli transpilation strategy to [PNNL NWQEC](https://github.com/pnnl/nwqec).

It is written so that another coding agent can begin implementation without reconstructing the design from the literature or rediscovering NWQEC's current structure. It specifies:

- the intended quantum transformation;
- the assumptions under which it is exact;
- the required changes to NWQEC's intermediate representation, passes, Python API, output, analysis, and tests;
- a safe staged implementation order;
- verification criteria that detect relative-phase and dynamic-control errors;
- expected gate counts and acceptance criteria.

Repository state referenced by this plan: the public `main` branch inspected on August 24, 2026. Before editing, the implementing agent should synchronize with the intended development branch and re-check whether the named files or APIs have changed.

## Executive summary

NWQEC currently lowers `CCX` to the standard exact circuit with seven \(T/T^\dagger\) gates and six `CX` gates. It also already supports `RCCX`, whose current decomposition uses four \(T/T^\dagger\) gates and three `CX` gates. The four-\(T\) `RCCX` is not, by itself, an exact Toffoli on arbitrary states because it introduces relative phases.

The proposed implementation uses one clean ancilla, one mid-circuit measurement, and a classically conditioned Clifford correction to turn the existing four-\(T\) `RCCX` into an exact deterministic `CCX` channel. The lowering is

```text
RCCX(control_0, control_1, ancilla)
SDG(ancilla)
CX(ancilla, target)
H(ancilla)
MEASURE(ancilla -> result_bit)
if result_bit == 1:
    CZ(control_0, control_1)
RESET(ancilla)  # only when the ancilla will be reused
```

For NWQEC's current `RCCX` phase convention, `SDG`, not `S`, is the required phase repair. This must be protected by a regression test.

The implementation should be introduced as a new semantic lowering pass, for example `CCXLoweringPass`, before the existing generic `DecomposePass`. The first release should support the dynamic four-\(T\) strategy in `to_clifford_t`, make it explicitly opt-in, and reject unsupported combinations such as dynamic lowering followed by the current unitary-only PBC/Tfuse pipeline. Later work can add dynamic-region-aware PBC optimization, ancilla-pool scheduling, temporary logical-AND recognition, and multi-controlled Toffoli lowering.

## Goals

The implementation must provide all of the following:

1. An exact logical `CCX`, not merely a relative-phase Toffoli.
2. Four total `T` or `TDG` operations per independently lowered `CCX`.
3. One clean ancilla per simultaneously active gadget, or one reusable clean ancilla when serialization is acceptable.
4. Mid-circuit measurement and a measurement-conditioned `CZ` correction.
5. Correct classical dependencies in scheduling and depth analysis.
6. A user-visible strategy selection API that preserves current behavior by default.
7. Verification of the induced quantum channel, including both measurement branches and superposition inputs.
8. Safe behavior when a downstream pass does not support dynamic circuits.

## Non-goals for the first implementation

The minimum viable implementation does not need to provide:

- automatic recognition of all compute/uncompute Toffoli pairs;
- measurement-aware Tfuse across dynamic boundaries;
- optimal ancilla scheduling for arbitrary parallel circuits;
- arbitrary multi-controlled `C^nX` synthesis;
- physical surface-code scheduling of measurement and feed-forward;
- a proof that four \(T\) gates is globally optimal under every possible catalytic resource model;
- automatic replacement of `CCX` by `RCCX` when relative phases happen to be harmless.

Those are follow-on features described later in this document.

## Background and exactness assumptions

### Static versus dynamic Toffoli costs

An exact static Clifford+T decomposition of a three-qubit Toffoli uses seven \(T/T^\dagger\) gates. Four-\(T\) unitary circuits commonly called `RCCX`, relative-phase Toffoli, Toffoli*, or logical AND do not generally implement the same unitary as `CCX` on arbitrary superpositions.

Cody Jones showed that an exact deterministic Toffoli can be implemented with four \(T\) gates if a clean ancilla, measurement, and classically controlled Clifford correction are allowed. The relevant reference is:

- Cody Jones, ["Low-overhead constructions for the fault-tolerant Toffoli gate," Physical Review A 87, 022328 (2013)](https://doi.org/10.1103/PhysRevA.87.022328).

Related references are:

- Peter Selinger, ["Quantum circuits of T-depth one," Physical Review A 87, 042302 (2013)](https://doi.org/10.1103/PhysRevA.87.042302).
- Dmitri Maslov, ["Advantages of using relative-phase Toffoli gates with an application to multiple control Toffoli optimization," Physical Review A 93, 022311 (2016)](https://doi.org/10.1103/PhysRevA.93.022311).
- Craig Gidney, ["Halving the cost of quantum addition," Quantum 2, 74 (2018)](https://doi.org/10.22331/q-2018-06-18-74).

### NWQEC's current `RCCX` convention

At the inspected repository revision, `include/nwqec/passes/decompose_pass.hpp` lowers `RCCX(c0,c1,t)` chronologically as

```text
H(t)
T(t)
CX(c1,t)
TDG(t)
CX(c0,t)
T(t)
CX(c1,t)
TDG(t)
H(t)
```

On computational-basis inputs, its nontrivial phase behavior includes

```text
|100> -> -|101> only when the target starts in 1 in the corresponding case
|110> ->  i|111>
|111> -> -i|110>
```

The only property needed by the ancilla-assisted construction is its action when the `RCCX` target is a clean ancilla initialized to \(|0\rangle\):

\[
|c_0c_1\rangle |0\rangle_a
\longmapsto
i^{c_0c_1}|c_0c_1\rangle |c_0c_1\rangle_a.
\]

Applying \(S^\dagger\) to the ancilla removes the phase because \(S^\dagger|1\rangle=-i|1\rangle\):

\[
S_a^\dagger\operatorname{RCCX}(c_0,c_1;a)
|c_0c_1\rangle|0\rangle_a
=
|c_0c_1\rangle|c_0c_1\rangle_a.
\]

Do not replace `SDG` with `S` by copying a circuit written under a different gate or phase convention. Verify the actual NWQEC sequence as part of the tests.

### Correctness of measurement-assisted uncomputation

After the corrected AND has been stored in the ancilla, apply `CX(ancilla,target)`. The logical data have now undergone the correct Toffoli permutation, but the ancilla remains entangled with the controls. Measuring the ancilla in the \(X\) basis is implemented by `H` followed by a computational-basis measurement.

Let the measurement outcome be \(m\in\{0,1\}\). Before feed-forward correction, the data acquire the phase

\[
(-1)^{m c_0c_1},
\]

which is exactly \(CZ(c_0,c_1)^m\). Applying a conditional `CZ` removes it. Therefore the corrected branch operator is

\[
K_m^{\mathrm{corrected}}
=
\frac{1}{\sqrt{2}}U_{\mathrm{CCX}}
\]

for both \(m=0\) and \(m=1\). Consequently,

\[
\sum_{m=0}^1
K_m^{\mathrm{corrected}}\rho
\left(K_m^{\mathrm{corrected}}\right)^\dagger
=
U_{\mathrm{CCX}}\rho U_{\mathrm{CCX}}^\dagger.
\]

Each measurement outcome has probability \(1/2\), independent of the input state. The implementation is deterministic after feed-forward; it is not a postselected or approximate gate.

## Current NWQEC structure relevant to this work

The implementing agent should inspect at least the following files before editing:

- `include/nwqec/core/operation.hpp`
  - Defines `Operation::Type`, including `CCX`, `RCCX`, `MEASURE`, `RESET`, and `CZ`.
  - Stores qubits, parameters, and classical measurement bits.
  - Does not currently store an execution condition on an operation.
  - Prints operations in OpenQASM-like syntax.
- `include/nwqec/core/circuit.hpp`
  - Owns flattened operations, qubit count, and classical-bit count.
  - Contains `is_clifford_t_operation`, gate counting, printing, depth, and duration methods.
  - Current depth analysis follows quantum operands but not classical feed-forward dependencies.
- `include/nwqec/passes/decompose_pass.hpp`
  - Contains the current seven-\(T\) `CCX` decomposition.
  - Contains the four-\(T\), three-`CX` `RCCX` decomposition used by this proposal.
  - Rebuilds a new circuit with the original number of qubits and bits.
- `include/nwqec/core/transpiler_passes.hpp`
  - Defines pass types and standard pass sequences.
- `include/nwqec/core/transpiler.hpp`
  - Runs configured pass sequences.
- `python/nwqec/_core.cpp`
  - Defines pybind11 bindings and top-level transforms such as `to_clifford_t`.
- `docs/python_api.md`
  - Documents the current public Python API.
- `tools/`
  - Contains the command-line entry points and will need equivalent strategy flags if the feature is exposed through the CLI.
- `tests/`
  - Add unit, integration, resource-count, serialization, and randomized channel tests here, following the repository's existing test conventions.

The current public API uses `keep_ccx: bool`. Backward compatibility should be maintained while introducing a more expressive strategy option.

## Proposed compiler architecture

### Add a semantic `CCXLoweringPass`

Create a new pass, suggested path:

```text
include/nwqec/passes/ccx_lowering_pass.hpp
```

This pass should decide how a high-level `CCX` is realized. It should run before the existing generic `DecomposePass`.

Suggested pipeline:

```text
Parse or construct Circuit
    -> CCXLoweringPass
    -> DecomposePass
    -> RZ synthesis and cleanup
    -> supported optimization passes
    -> backend-specific dynamic lowering or output
```

Do not put ancilla allocation and dynamic-control policy directly inside `DecomposePass`. `DecomposePass` should remain responsible for local gate-to-basis translation. The new pass is responsible for semantic strategy and resources.

### Strategy and options types

Add an enum and options object similar to:

```cpp
enum class CCXStrategy {
    Preserve,
    Standard7T,
    Jones4TDynamic,
    RelativePhase4T,
    Auto
};

struct CCXLoweringOptions {
    CCXStrategy strategy = CCXStrategy::Standard7T;
    size_t max_clean_ancillas = 0;
    bool reuse_ancillas = true;
    bool emit_reset = true;
    bool allow_measurement = false;
    bool allow_classical_feedforward = false;
};
```

Required meanings:

- `Preserve`: leave `CCX` unchanged.
- `Standard7T`: current exact unitary lowering; no ancilla or measurement.
- `Jones4TDynamic`: exact four-\(T\) lowering described here.
- `RelativePhase4T`: explicitly request a relative-phase replacement. This must never be labeled exact `CCX` lowering.
- `Auto`: select a strategy using backend capability and resource constraints. In the first release, it may simply select `Jones4TDynamic` only when dynamic execution is explicitly enabled and at least one clean ancilla is available; otherwise select `Standard7T`.

`Jones4TDynamic` must fail with a clear diagnostic when measurement, feed-forward, or a clean ancilla is unavailable. Do not silently fall back unless an explicit option requests fallback behavior.

## Intermediate-representation changes

### Add a generic classical condition to `Operation`

The key missing IR feature is a condition attached to an otherwise ordinary operation. Add a structure similar to:

```cpp
struct ClassicalCondition {
    std::vector<size_t> bits;
    uint64_t value = 0;

    bool operator==(const ClassicalCondition&) const = default;
};
```

Extend `Operation` with:

```cpp
std::optional<ClassicalCondition> condition;
```

Add accessors such as:

```cpp
bool is_conditional() const;
const std::optional<ClassicalCondition>& get_condition() const;
Operation with_condition(ClassicalCondition condition) const;
```

Implementation requirements:

- Existing constructors must remain source compatible where practical.
- Copying an `Operation` must preserve its condition.
- Any helper that manually reconstructs an operation must copy the condition, dagger flag, `x_rotation`, Pauli metadata, parameters, qubits, and classical bits.
- `Circuit::expand_gate` currently reconstructs operations from selected fields; audit it and other similar sites so new metadata is not dropped.
- Equality, hashing, debugging, Python exposure, and serialization should include the condition where applicable.

Do not introduce a special `CONDITIONAL_CZ` gate unless the existing architecture makes generic conditioning impossible. Generic conditioning will be useful for teleportation, active reset, magic-state injection, and future dynamic protocols.

### Classify dynamic operations

Add helpers such as:

```cpp
bool Operation::is_measurement() const;
bool Operation::is_reset() const;
bool Operation::is_conditional() const;
bool Operation::is_dynamic() const;

bool Circuit::has_measurement() const;
bool Circuit::has_feedforward() const;
bool Circuit::is_dynamic() const;
```

These should be used by pipeline validation rather than repeatedly scanning for individual operation types in unrelated passes.

### Preserve classical dependencies

A conditional operation reads its condition bits. A measurement writes its destination bits. Represent that relationship in dependency analysis.

At minimum, dependency construction must add

```text
MEASURE(ancilla -> result_bit)
    -> CONDITIONAL CZ(control_0, control_1) conditioned on result_bit
```

No scheduler or optimization pass may move the conditional correction before its measurement.

## Ancilla and classical-bit management

### Minimal allocator

Add an allocator owned by the lowering pass or a reusable circuit-resource utility:

```cpp
class AncillaManager {
public:
    size_t acquire_clean_qubit();
    void release_clean_qubit(size_t qubit);

    size_t acquire_classical_bit();
    void release_classical_bit(size_t bit);
};
```

The allocator must know:

- which qubits were part of the input circuit;
- which appended qubits are compiler-created clean ancillas;
- whether a released ancilla has been reset and is safe to reuse;
- which classical bits store live measurement outcomes;
- the configured maximum ancilla count.

### Resource semantics

A newly allocated clean ancilla is assumed to begin in \(|0\rangle\). If the same physical ancilla is reused after measurement, insert `RESET` unless the backend contract guarantees a fresh clean replacement.

Never treat an arbitrary unused or dirty data qubit as clean.

For the first version, one reusable ancilla is acceptable, but document that it serializes all lowered Toffolis through a shared resource. A later pool-aware scheduler should support:

- one shared ancilla, minimizing space;
- one ancilla per potentially parallel gadget, maximizing available parallelism;
- a fixed pool of \(k\) ancillas, trading space for depth.

### Circuit register handling

`DecomposePass` currently creates flattened `q` and `c` registers sized from the source circuit. The new pass must increase both counts safely.

Recommended first-version behavior:

1. Determine the required ancilla pool size and number of reusable result bits before rebuilding the output circuit.
2. Create the output quantum and classical register sizes once.
3. Assign appended global indices to compiler-created resources.
4. Preserve original qubit and bit indices.
5. Ensure printed register sizes, register metadata, `get_num_qubits`, and `get_num_bits` agree.

Do not rely only on `add_operation` implicitly increasing the maximum index if that leaves register-size maps inconsistent.

## Exact lowering algorithm

### Input validation

For every `CCX` operation selected for `Jones4TDynamic`:

1. Verify it has exactly three distinct qubit operands.
2. Interpret operands using NWQEC's current ordering: control 0, control 1, target.
3. Acquire a clean ancilla distinct from all data operands.
4. Acquire a classical result bit.
5. Verify that dynamic execution and classical feed-forward are enabled.

### Emitted high-level operations

Emit, in order:

```cpp
Operation(RCCX,   {control_0, control_1, ancilla});
Operation(SDG,    {ancilla});
Operation(CX,     {ancilla, target});
Operation(H,      {ancilla});
Operation(MEASURE,{ancilla}, /* parameters */ {}, {result_bit});

Operation correction(CZ, {control_0, control_1});
correction = correction.with_condition({{result_bit}, 1});

if (reuse_ancillas && emit_reset) {
    Operation(RESET, {ancilla});
}
```

Check the actual `Operation` constructor argument order before implementing; the pseudocode conveys semantics, not a guaranteed compilable signature.

### Resource counts before backend lowering

After generic decomposition of `RCCX`, but before physically decomposing the conditional `CZ`, each lowered `CCX` should contain:

- four `T` or `TDG` gates total;
- three `CX` gates inside `RCCX`;
- one additional `CX(ancilla,target)`;
- two `H` gates inside `RCCX`;
- one additional `H` for \(X\)-basis measurement;
- one `SDG` phase correction;
- one measurement;
- one conditional `CZ`;
- zero or one reset depending on reuse policy.

Thus the base dynamic gadget has:

\[
N_T=4,\qquad N_{CX}=4,
\]

where the conditional `CZ` remains a logical Clifford correction. If it is decomposed physically as `H-CX-H`, the physical count becomes five `CX` gates and five `H` gates in total.

Gate-count tests must state clearly whether conditional corrections are retained, physically decomposed, or tracked in a Clifford frame.

### Conditional correction lowering

Support two backend choices:

1. Physical dynamic correction:
   - retain a conditional `CZ`, or decompose it into conditional basis gates;
   - all decomposed components must inherit the same classical condition.
2. Conditional Clifford-frame update:
   - record the `CZ` byproduct in a classically conditioned Clifford frame;
   - do not count it as an executed quantum gate;
   - preserve it through later logical operations using correct frame propagation.

The first implementation should prefer the physical dynamic correction unless NWQEC already has a branch-conditioned Clifford-frame abstraction. Do not silently discard the `CZ` merely because it has zero \(T\) cost.

## Changes to the generic decomposition pass

The existing `DecomposePass` should keep its current `RCCX` decomposition. Make only the changes required to preserve dynamic metadata:

- If an input operation is conditional and is decomposed into several operations, every emitted operation must inherit the same condition.
- Measurements and resets should be retained as supported basis operations.
- A conditional `CZ` should either remain intact for a dynamic backend or lower to conditional `H`, `CX`, `H` operations.
- Unsupported dynamic lowering must produce a diagnostic rather than an unconditional circuit.

Retain the standard `CCX` decomposition for `Standard7T`. Move the choice between preserving, standard lowering, dynamic lowering, and relative-phase lowering into `CCXLoweringPass`.

## Pass pipeline and compatibility policy

### Clifford+T pipeline

The minimum viable feature should work in `to_clifford_t`:

```text
CCXLoweringPass(strategy)
DecomposePass
RZ synthesis if needed
safe cleanup passes
```

The output is Clifford+T plus measurement, reset, and classical feed-forward.

### PBC and Tfuse

The current PBC/Tfuse flow should be treated as unitary-only unless its implementation explicitly proves otherwise.

For the first version:

- reject `Jones4TDynamic` combined with `to_pbc`, `fuse_t`, or any pass that cannot preserve measurements and conditional branches;
- emit a precise error such as: `Jones4TDynamic introduces measurement and feed-forward; the current PBC pipeline supports only unitary regions`;
- allow `Standard7T` to use the existing PBC flow;
- allow a future region-based pipeline to optimize unitary segments separated by dynamic barriers.

Do not lower every `CCX` dynamically before global phase-polynomial optimization without considering the effect. Measurement barriers can prevent optimizations that would combine phases across neighboring gates. A future `Auto` cost model should compare at least:

```text
standard 7-T lowering followed by global optimization
```

against

```text
local 4-T dynamic lowering with measurement barriers
```

### Optimization barriers

Until a pass is explicitly dynamic-safe, treat the following as barriers:

- `MEASURE`;
- `RESET`;
- any conditional operation;
- any branch-conditioned frame update.

A pass may optimize within each unitary region but must not commute an operation across a barrier based only on qubit commutation.

## Python API and backward compatibility

### Proposed API

Expose a string or enum-like strategy through pybind11:

```python
nwqec.to_clifford_t(
    circuit,
    ccx_strategy="jones-4t",
    max_clean_ancillas=1,
    reuse_ancillas=True,
    allow_dynamic=True,
)
```

Supported values should be documented as:

```text
"preserve"
"standard-7t"
"jones-4t"
"relative-phase-4t"
"auto"
```

Use the same options in `get_clifford_t_counts` and, when appropriate, the CLI.

### Backward compatibility

Keep `keep_ccx` temporarily and map it as follows:

```text
keep_ccx=True  -> ccx_strategy="preserve"
keep_ccx=False -> ccx_strategy="standard-7t"
```

If both the old and new arguments are supplied inconsistently, raise an error instead of guessing.

The default must remain `standard-7t` for the first release so existing callers do not unexpectedly receive measurements, extra qubits, or conditional control.

### Counts and statistics

The current count API returns a gate-count dictionary. Dynamic compilation also needs resource metadata. Either extend the returned data compatibly or add a new result type such as:

```python
result.gate_counts
result.num_data_qubits
result.num_clean_ancillas
result.num_classical_bits_added
result.num_measurements
result.num_resets
result.num_conditional_operations
result.t_count
result.t_depth
result.adaptive_depth
```

Do not report only `T=4` while hiding that the implementation requires an ancilla, measurement, and feed-forward.

## QASM input and output

### OpenQASM 2 limitation

OpenQASM 2 conditions are expressed against an entire classical register, for example:

```qasm
if (c == 1) cz q[0], q[1];
```

NWQEC currently flattens classical storage into `creg c[n]`. A condition intended to read one bit cannot generally be printed as `if (c == 1)` when other bits coexist in `c`.

Do not emit ambiguous or semantically incorrect OpenQASM 2.

### Recommended output solution

Prefer adding an OpenQASM 3 writer for dynamic circuits:

```qasm
bit result;
result = measure ancilla;
if (result == 1) {
    cz control_0, control_1;
}
```

If OpenQASM 2 output must be supported, preserve or create separate one-bit classical registers for independent conditions and ensure the register comparison is exact. This likely requires preserving register identity instead of flattening every bit into a single printed register.

The internal IR should still use individual bit indices, independent of output language.

## Scheduling, depth, and duration

### Classical dependency edges

Extend depth and DAG construction so a conditional operation depends on the last writer of every condition bit.

Track at least:

```text
last quantum operation on each qubit
last writer of each classical bit
classical bits read by each conditional operation
```

For a conditional operation, its earliest start is constrained by both its qubit dependencies and its condition-bit dependencies.

### Latency model

Add configurable durations for:

- measurement;
- classical decoding or feed-forward;
- conditional Clifford execution or frame update;
- reset.

At minimum, provide separate metrics:

```text
quantum gate depth
T depth
measurement depth
adaptive depth
estimated duration
```

The existing `Circuit::depth()` may remain as a simple gate-layer metric for backward compatibility, but a dynamic-aware depth function must be added or the old function must be corrected and documented.

### Ancilla reuse and parallelism

One shared ancilla creates a dependency chain through all dynamic Toffolis and can significantly increase depth. Make this visible in statistics.

For a later pool-aware implementation, assign gadgets to ancillas using resource availability rather than simply releasing and immediately reusing the same ancilla while scanning a linear operation list. A practical heuristic is to assign each new gadget to the ancilla with the smallest scheduled availability time, subject to `max_clean_ancillas`.

## Testing and formal verification plan

Computational-basis truth-table tests are necessary but not sufficient. A relative-phase Toffoli passes every classical truth-table test while still being the wrong quantum gate.

### Test level 1: current `RCCX` phase convention

Construct the matrix of NWQEC's exact emitted `RCCX` sequence and verify all eight computational-basis transitions including phases.

At minimum verify:

\[
\operatorname{RCCX}|110\rangle=i|111\rangle,
\qquad
\operatorname{RCCX}|111\rangle=-i|110\rangle.
\]

Then verify the clean-target identity:

\[
S_a^\dagger\operatorname{RCCX}(c_0,c_1;a)
|c_0c_1\rangle|0\rangle_a
=
|c_0c_1\rangle|c_0c_1\rangle_a
\]

for all four control inputs.

This test is the guard against changing `SDG` to `S` or changing the `RCCX` decomposition without updating the dynamic lowering.

### Test level 2: basis-state logical behavior

For all eight data inputs \(|c_0c_1t\rangle\):

1. initialize the ancilla to \(|0\rangle\);
2. simulate the dynamic gadget;
3. evaluate both measurement outcomes;
4. apply the appropriate correction;
5. verify that the data output is \(|c_0,c_1,t\oplus c_0c_1\rangle\);
6. verify that each branch has probability \(1/2\).

### Test level 3: branchwise Kraus-operator equality

Build the data-space Kraus operator for each corrected measurement branch and verify

\[
K_0=\frac{1}{\sqrt{2}}U_{\mathrm{CCX}},
\qquad
K_1=\frac{1}{\sqrt{2}}U_{\mathrm{CCX}}.
\]

Compare complex matrices within a strict numerical tolerance such as \(10^{-12}\). If the simulator introduces a harmless global phase per branch, normalize it before comparison and separately verify equal branch probabilities.

This is the primary exactness test.

### Test level 4: channel equality on random states

Generate random normalized complex three-qubit states \(|\psi\rangle\), not only real-valued states. Verify

\[
\sum_m K_m|\psi\rangle\!\langle\psi|K_m^\dagger
=
U_{\mathrm{CCX}}|\psi\rangle\!\langle\psi|U_{\mathrm{CCX}}^\dagger.
\]

Use a fixed random seed and multiple trials. Include deliberately phase-sensitive states such as

\[
\frac{|000\rangle+|110\rangle}{\sqrt{2}},
\qquad
\frac{|101\rangle+i|110\rangle}{\sqrt{2}}.
\]

The first of these will expose many relative-phase mistakes immediately.

### Test level 5: complete-channel or Choi-state equality

To test the action on data entangled with an external reference, prepare a maximally entangled state between the three data qubits and three reference qubits. Apply the compiled channel to the data half and compare with the ideal `CCX` channel.

Equivalent alternatives are:

- compare Choi matrices;
- compare process matrices;
- verify the diamond distance numerically for this small channel, if a trusted library is available.

Choi equality is strongly recommended because it verifies the complete quantum channel, not just its action on selected unentangled states.

### Test level 6: resource counts

For one `CCX` compiled with `Jones4TDynamic`, assert:

```text
T + TDG count                    = 4
base CX count                    = 4
SDG count                        = 1
measurement count                = 1
conditional CZ count             = 1
clean ancillas added              = 1
classical result bits added       = 1, or the documented reusable count
reset count                       = 1 when reuse/reset is enabled
```

If the conditional `CZ` is decomposed physically, assert the documented post-lowering counts separately. Do not mix logical and physical correction counts in one expected value.

Also verify the unchanged strategies:

```text
Standard7T: exactly 7 T/TDG, no added ancilla, no measurement
Preserve: one CCX remains
RelativePhase4T: exactly 4 T/TDG and is explicitly marked non-equivalent to exact CCX
```

### Test level 7: multiple gates and ancilla reuse

Test circuits containing:

- two sequential `CCX` gates with the same controls;
- two `CCX` gates on disjoint data qubits;
- `CCX` gates separated by Clifford gates;
- `CCX` gates separated by non-Clifford gates;
- a mix of preserved, standard, and dynamic strategies where the API permits it;
- pre-existing measurements and classical bits.

Verify:

- the ancilla is reset before reuse;
- result bits are not overwritten before the associated correction consumes them;
- one shared ancilla introduces the expected dependency;
- an ancilla pool does not introduce false data dependencies;
- the overall logical channel equals the ideal circuit.

### Test level 8: serialization and parser/writer round trips

For dynamic output:

1. compile a circuit with one `CCX`;
2. serialize it;
3. parse it again where supported;
4. verify the measurement destination and correction condition survive;
5. compare the round-tripped logical channel;
6. ensure individual-bit conditions are not incorrectly widened to an entire multi-bit register.

If OpenQASM 2 cannot represent the exact internal condition safely, the writer must reject the output with a clear message rather than generating incorrect QASM.

### Test level 9: unsupported-pipeline failures

Add negative tests that verify informative failure for:

- `Jones4TDynamic` with `allow_dynamic=False`;
- no available clean ancilla;
- feed-forward-disabled backend;
- current unitary-only PBC/Tfuse after dynamic lowering;
- QASM 2 output when the condition cannot be represented exactly;
- malformed `CCX` operands;
- inconsistent `keep_ccx` and `ccx_strategy` arguments.

### Test level 10: regression and performance

Run the complete existing test suite. Confirm that default compilation still uses the current standard seven-\(T\) circuit and produces unchanged results.

For benchmark circuits with many Toffolis, report:

- old and new \(T\) counts;
- old and new data-plus-ancilla qubit counts;
- measurement and reset counts;
- ordinary depth and adaptive depth;
- output size;
- transpilation time.

Do not claim a runtime improvement solely from the reduction \(7\to4\) without including measurement, feed-forward, reset, and serialization-induced loss of parallelism.

## Suggested test implementation utilities

Create a small independent simulator or matrix utility in the test suite rather than validating the pass with the same algebra used to implement it.

Useful helpers include:

```python
def operation_matrix(operation, num_qubits): ...
def unitary_of_static_region(operations, num_qubits): ...
def branch_kraus_operators(dynamic_circuit, measured_bit): ...
def ideal_ccx_matrix(): ...
def choi_matrix(kraus_ops): ...
def equal_up_to_global_phase(a, b, atol=1e-12): ...
```

For the dynamic gadget, an even safer test is to construct branch Kraus operators directly from projectors on the measured ancilla:

\[
K_m
=
\langle m|_a
H_a
CX_{a,t}
S_a^\dagger
RCCX_{c_0,c_1;a}
|0\rangle_a,
\]

followed by \(CZ_{c_0,c_1}^m\).

This test should not depend on NWQEC's QASM writer or its classical-control simulator.

## Documentation requirements

Update at least:

- `README.md` quick-start example;
- `docs/python_api.md`;
- CLI documentation if a CLI flag is added;
- release notes or changelog;
- source-level comments in `CCXLoweringPass` citing Jones 2013;
- a short explanation distinguishing exact `CCX` from `RCCX`;
- a warning that the dynamic strategy changes qubit, classical-bit, measurement, and latency requirements.

Suggested example:

```python
import nwqec

circuit = nwqec.load_qasm("input.qasm")

compiled = nwqec.to_clifford_t(
    circuit,
    ccx_strategy="jones-4t",
    max_clean_ancillas=1,
    reuse_ancillas=True,
    allow_dynamic=True,
)

print(compiled.stats())
```

The documentation must state that `jones-4t` produces an exact dynamic circuit, while `relative-phase-4t` is not a drop-in exact replacement for arbitrary `CCX` gates.

## Phased implementation tasks

### Phase 0: baseline and convention capture

- Check out the intended branch and run the full existing test suite.
- Record the baseline result.
- Add a test that captures the current seven-\(T\) `CCX` count.
- Add a matrix test that captures the current `RCCX` phase convention.
- Confirm operand order for `CCX`, `RCCX`, `CX`, measurement bits, and `CZ`.

Exit criteria:

- Existing tests pass.
- The `RCCX|110\rangle=i|111\rangle` convention is verified independently.

### Phase 1: classical-condition IR

- Add `ClassicalCondition` and optional condition metadata to `Operation`.
- Update constructors, accessors, copies, reconstruction sites, and Python bindings.
- Add dynamic-circuit query helpers.
- Add printing/debug representation without yet promising full QASM round-trip support.

Exit criteria:

- An operation can be conditioned on an individual bit.
- The condition survives circuit copies and pass-through passes.
- Existing unconditional circuits remain unchanged.

### Phase 2: exact lowering with one newly allocated ancilla per `CCX`

- Add `CCXStrategy` and `CCXLoweringPass`.
- Initially allocate a distinct ancilla and result bit for each `CCX`; this is easier to verify and does not create accidental cross-gadget dependencies.
- Emit the exact sequence with `SDG`.
- Run existing `RCCX` decomposition afterward.
- Add branchwise Kraus and resource-count tests.

Exit criteria:

- Every corrected branch equals \(U_{CCX}/\sqrt{2}\).
- Each gadget contains four \(T/T^\dagger\) gates.
- No relative-phase-only behavior remains on the data.

### Phase 3: reusable ancilla and result-bit pool

- Add the resource manager.
- Support a pool size set by `max_clean_ancillas`.
- Insert reset before reuse where required.
- Prevent premature classical-bit reuse.
- Add scheduling and reuse tests.

Exit criteria:

- A pool of one works correctly for multiple Toffolis.
- A larger pool can preserve some independent-gate parallelism.
- Resource statistics report the pool accurately.

### Phase 4: Python API, counts, CLI, and documentation

- Expose options through pybind11.
- Preserve `keep_ccx` compatibility.
- Extend count and statistics APIs.
- Add CLI flags if the CLI supports Clifford+T conversion options.
- Update documentation and examples.

Exit criteria:

- Python users can explicitly select all strategies.
- Default behavior is unchanged.
- Conflicting or unsupported options fail clearly.

### Phase 5: dynamic-aware scheduling and output

- Add classical dependency edges and adaptive latency.
- Implement correct dynamic serialization, preferably OpenQASM 3.
- Add parser/writer round-trip tests where supported.
- Decide whether conditional `CZ` remains physical or becomes a Clifford-frame update.

Exit criteria:

- The correction cannot be scheduled before measurement.
- Output preserves the exact individual-bit condition.
- Reported adaptive duration includes measurement and feed-forward.

### Phase 6: safe optimization integration

- Mark dynamic boundaries explicitly.
- Reject unsupported PBC/Tfuse use in the meantime.
- Optionally partition circuits into unitary regions and optimize each region independently.
- Add equivalence tests around every optimization path.

Exit criteria:

- No existing pass drops, duplicates, or moves a measurement or conditional correction incorrectly.
- Dynamic compilation either succeeds correctly or fails explicitly.

## Follow-on optimization: temporary logical AND

After the isolated exact four-\(T\) Toffoli works, the next high-value feature is a `TemporaryANDPass` based on Gidney's measurement-based uncomputation.

Many arithmetic circuits use a pattern equivalent to

```text
CCX(a, b, work)       # compute work = a AND b
... use work ...
CCX(a, b, work)       # uncompute work
```

Compiling both gates independently with Jones gadgets costs eight \(T\) gates. A temporary logical AND can compute with four \(T\) gates and erase with measurement and zero additional \(T\) gates, reducing the pair to four total.

Pattern recognition must initially be conservative. Required proof obligations include:

- `work` is known to start in \(|0\rangle\);
- the first operation semantically computes the AND;
- the intermediate uses of `work` are compatible with measurement-based erasure, commonly control-only uses;
- no intervening operation requires preservation of an untracked relative phase;
- the final operation solely uncomputes the work bit;
- measurement and correction dependencies are representable.

A safer first interface may introduce explicit semantic operations such as `AND_COMPUTE` and `AND_UNCOMPUTE`, then add automatic recognition later.

## Follow-on optimization: multi-controlled Toffoli

Once clean and conditionally clean ancilla management exists, extend the strategy system to `C^nX` gates. Relevant recent references include:

- Junhong Nie, Wei Zi, and Xiaoming Sun, ["Quantum circuit for multi-qubit Toffoli gate with optimal resource" (2024)](https://arxiv.org/abs/2402.05053).
- Suman Dutta et al., ["On Exact Space-Depth Trade-Offs in Multi-Controlled Toffoli Decomposition" (2025)](https://arxiv.org/abs/2502.01433).
- Abhoy Kole et al., ["Adaptive Clifford+T Decomposition of Large Toffoli Gates with One Clean Ancilla" (2026)](https://arxiv.org/abs/2605.18169).
- Abhoy Kole et al., ["Measurement-Driven Adaptive Low-Overhead Implementation of Multi-Controlled Toffoli Gates" (2026)](https://arxiv.org/abs/2605.18159).

This should be a separate follow-on because its optimal choice depends strongly on ancilla count, conditional cleanliness, connectivity, T-depth, CX overhead, and measurement latency.

## Risks and mitigations

### Risk: confusing `RCCX` with exact `CCX`

Mitigation:

- keep distinct strategy names;
- require an explicit relative-phase strategy;
- verify superposition states and Kraus operators;
- document exactness prominently.

### Risk: wrong `S` versus `SDG` phase repair

Mitigation:

- derive from NWQEC's emitted `RCCX`, not a diagram from another convention;
- lock the phase convention with a matrix regression test;
- test both measurement branches.

### Risk: dropping the conditional correction

Mitigation:

- represent conditions in the IR;
- classify conditional operations as dynamic barriers;
- add Choi-channel equality tests;
- fail unsupported backends explicitly.

### Risk: invalid QASM 2 condition

Mitigation:

- prefer OpenQASM 3 for dynamic output;
- preserve one-bit registers if QASM 2 is required;
- test serialization semantically;
- reject unrepresentable conditions.

### Risk: apparent T-count win but worse runtime

Mitigation:

- report ancillas, measurements, adaptive depth, resets, and feed-forward latency;
- provide standard and dynamic strategies;
- make `Auto` backend- and cost-model-aware in later releases.

### Risk: optimization across a measurement boundary

Mitigation:

- insert explicit barriers in dependency analysis;
- reject dynamic input in unverified passes;
- optimize only within proven unitary regions.

### Risk: shared ancilla destroys parallelism

Mitigation:

- expose pool size;
- report depth impact;
- later assign ancillas by scheduled availability rather than source-list order.

## Definition of done for the first production-ready release

The feature is complete only when all of the following are true:

- `to_clifford_t` accepts an explicit exact four-\(T\) dynamic CCX strategy.
- The old default behavior remains the standard seven-\(T\) decomposition.
- The implementation uses the existing `RCCX` plus the correct `SDG` repair.
- One clean ancilla and measurement outcome are represented explicitly.
- The `CZ` correction is conditioned on the correct individual bit.
- Both measurement branches produce the ideal Toffoli channel.
- Four total \(T/T^\dagger\) gates are emitted per independently lowered `CCX`.
- Ancilla and classical-bit counts are reported.
- Reuse includes reset and correct dependencies.
- Dynamic-unsafe passes reject the circuit clearly.
- Serialization either preserves exact semantics or explicitly rejects unsupported output.
- Existing NWQEC tests pass without changing default outputs.
- New basis, phase, Kraus, Choi, resource, reuse, failure, and round-trip tests pass.
- User documentation explains the resource and latency tradeoffs.

## Recommended first coding task

The first coding agent should begin with Phase 0 and Phase 1 only:

1. run the existing test suite and record the baseline;
2. add an independent matrix test for the current `RCCX` phase convention;
3. add a `ClassicalCondition` field to `Operation`;
4. ensure it survives copying, reconstruction, Python binding, and printing;
5. add a unit test containing `MEASURE(a -> m)` followed by `CZ(c0,c1)` conditioned on `m`;
6. stop before implementing the lowering if condition metadata is lost anywhere.

Once that infrastructure is reliable, the four-\(T\) quantum lowering is small and can be verified branchwise before integration into larger pipelines.

