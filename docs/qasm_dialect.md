QASM Dialect Reference
======================

NWQEC reads and writes a dialect of OpenQASM 2.0. This document lists exactly which gates are recognized, which parts of the dialect go beyond OpenQASM 2, and which parts of the standard library are not supported.

Every statement in this document was checked against the parser and the exporter rather than inferred from the source, so the tables can be read as behaviour, not intent.

How NWQEC Reads a File
----------------------

Three things differ from a strict OpenQASM 2.0 reader and are worth knowing before the tables below make sense.

- **The `include` line is ignored.** NWQEC resolves gate names against a fixed builtin table compiled into the parser, not against `qelib1.inc`. `include "qelib1.inc";` is accepted and discarded, and so is `include "anything-else.inc";` — a nonexistent file raises no error. The consequence is that the supported gate set is exactly the builtin table plus whatever the file defines itself, regardless of what is included.
- **The `OPENQASM 2.0;` header is optional.** A file that begins directly with `qreg` parses.
- **User `gate` definitions are inlined at parse time.** A `gate` block is expanded into its body when applied, and the definition itself does not survive into the output. This is the escape hatch for everything in the "Not Supported" section below.

Builtin Gate Set
----------------

The 44 gate names below are recognized without any definition in the file. The `qelib1` column records whether the name is also in the standard `qelib1.inc` shipped with Qiskit, which is the de-facto OpenQASM 2 standard library.

### Single-qubit

| name | arguments | qelib1 |
|---|---|---|
| `x`, `y`, `z`, `h` | — | yes |
| `s`, `sdg`, `t`, `tdg` | — | yes |
| `sx`, `sxdg` | — | yes |
| `id` | — | yes |
| `rx`, `ry`, `rz` | 1 | yes |
| `p`, `u1` | 1 | yes |
| `u2` | 2 | yes |
| `u`, `u3` | 3 | yes |

### Two-qubit

| name | arguments | qelib1 |
|---|---|---|
| `cx`, `cy`, `cz`, `ch` | — | yes |
| `csx` | — | yes |
| `swap` | — | yes |
| **`cs`, `csdg`, `ct`, `ctdg`** | — | **no** |
| **`ecr`** | — | **no** |
| `crx`, `cry`, `crz` | 1 | yes |
| `cp`, `cu1` | 1 | yes |
| `cu3` | 3 | yes |
| `cu` | 4 | yes |
| `rxx`, `rzz` | 1 | yes |
| **`ryy`** | 1 | **no** |

### Three-qubit and multi-controlled

| name | arguments | qelib1 |
|---|---|---|
| `ccx`, `cswap`, `rccx` | — | yes |
| **`mcx`** | — | **no** |

### Non-gate statements

`measure`, `reset`, `barrier`, `qreg`, `creg`, `gate`, and the `if` statement described below.

Extensions Beyond OpenQASM 2
----------------------------

### 1. `mcx` — variable-arity multi-controlled X

```qasm
mcx q[0],q[1],q[2],q[3],q[4];   // C^4X: four controls, target last
```

The operand list is `{c_0, ..., c_{n-1}, target}`, so the **last** operand is the target and the arity is unbounded. This is the only builtin whose operand count is not fixed. OpenQASM 2 has no equivalent; the standard library instead provides the fixed-size `c3x` and `c4x`, neither of which NWQEC recognizes.

`mcx` is preserved through the pipeline by default. To lower it, select a strategy — see `docs/python_api.md` for `mcx_strategy` and `mcx_max_ancillas`, and `results/260829_gidney_mcx/` for the resource trade-offs between them.

### 2. Extra two-qubit gates

`cs`, `csdg`, `ct`, `ctdg`, `ecr` and `ryy` are recognized but are not in `qelib1.inc`. A file using them is **not portable** to a strict OpenQASM 2 reader unless it also carries its own `gate` definitions for them. Note in particular that `qelib1.inc` defines `rxx` and `rzz` but *not* `ryy`, so the `ryy` asymmetry is a property of the standard library rather than an oversight here.

### 3. Classically conditioned operations — OpenQASM 3 `if`

This is the one place where NWQEC deliberately borrows syntax from OpenQASM 3 inside an otherwise OpenQASM 2 body. It is the dialect the transpiler emits whenever a circuit contains feed-forward, which includes every measurement-based ancilla erasure produced by the CCX and MCX lowering passes.

Four forms are accepted on input:

```qasm
if (c == 1) x q[0];                 // whole register, single statement
if (c[0] == 1) x q[0];              // single bit, single statement
if (c == 3) { x q[0]; h q[1]; }     // whole register, block
if (c[1] == 0) z q[0];              // single bit, testing zero
```

Output is always normalized to the braced form, and a file containing any conditional carries a marker comment in its preamble:

```qasm
OPENQASM 2.0;
include "qelib1.inc";
// NWQEC dynamic dialect: conditional statements use OpenQASM 3 if-syntax.
```

The condition compares a classical register, or one bit of it, against an integer. Adjacent operations sharing an identical condition are grouped into a single `if` block on output rather than each getting its own.

Two constraints: `measure`, `reset` and `barrier` cannot be conditioned, and a condition may span at most 64 classical bits.

A complete example, as emitted by the Gidney AND erasure in the MCX and CCX lowering passes:

```qasm
h q[5];
measure q[5] -> c[1];
if (c[1] == 1) { h q[2]; cx q[4],q[2]; h q[2]; }
reset q[5];
```

### 4. The PBC dialect

Circuits converted with `to_pbc` are written in a Pauli-based dialect that shares OpenQASM 2's file structure but none of its gate vocabulary. It is not readable by any other OpenQASM tool.

| statement | meaning |
|---|---|
| `t_pauli <sign><pauli>;` | rotation by $\pi/4$ about the given Pauli operator |
| `s_pauli <sign><pauli>;` | rotation by $\pi/2$ |
| `z_pauli <sign><pauli>;` | rotation by $\pi$ |
| `m_pauli <sign><pauli>;` | terminal projective Pauli measurement |
| `m_pauli <sign><pauli> -> c[i];` | mid-circuit Pauli measurement writing bit `i` |

The Pauli string has one character per qubit from `{I, X, Y, Z}`, prefixed by `+` or `-`; a `-` denotes the negated operator, and therefore a negated rotation angle. The presence or absence of the `-> c[i]` destination is what distinguishes a mid-circuit measurement from a terminal one on readback.

```qasm
t_pauli +XZI;
m_pauli +XII -> c[0];
if (c[0] == 1) { z_pauli +ZXI; z_pauli +XZI; }
t_pauli +IIZ;
m_pauli +XII;
```

Conditionals compose with the PBC dialect exactly as above, which is how feed-forward survives the conversion.

Lowering to Clifford+T
----------------------

Parsing a gate and *lowering* it are separate questions. Every builtin below was transpiled with `to_clifford_t` and the resulting unitary compared against Qiskit's reference for the same gate; the $T$-counts are measured, not derived.

### Exact, no approximation

These lower to Clifford+T with no synthesis error. $T$-counts are for a single application.

| gates | $T$ |
|---|---|
| `x`, `y`, `z`, `h`, `s`, `sdg`, `sx`, `sxdg`, `id` | 0 |
| `cx`, `cy`, `cz`, `swap`, `ecr` | 0 |
| `t`, `tdg` | 1 |
| `ch` | 2 |
| `cs`, `csdg`, `csx` | 3 |
| `rccx` | 4 |
| `ccx`, `cswap` | 7 |
| `mcx` | strategy-dependent; $4n-4$ with `and-ladder` |

`mcx` is the one entry that does **not** lower by default: its default strategy is `preserve`, so an `mcx` survives `to_clifford_t` untouched unless `mcx_strategy` is set. This is deliberate, since the strategies trade $T$-count against ancilla count and the choice belongs to the caller.

### Approximate, via gridsynth RZ synthesis

These contain arbitrary-angle rotations, so they are synthesized to within the `rz_err` / `epsilon` budget rather than decomposed exactly. Counts below are at the default `per-gate` policy with $\varepsilon = 10^{-10}$ and are angle-dependent; treat them as representative, not fixed.

| gates | $T$, representative |
|---|---|
| `rx`, `ry`, `rz`, `p`, `u1`, `rxx`, `ryy`, `rzz` | ~100 |
| `u2` | ~198 |
| `crx`, `cry`, `crz` | ~200 |
| `u`, `u3`, `ct` | ~296 |
| `ctdg` | ~298 |
| `cp`, `cu1` | ~300 |
| `cu3` | ~496 |
| `cu` | ~498 |

Spot-checked against Qiskit: `rz`, `rx`, `ryy`, `rxx` agree to $7 \times 10^{-11}$; `u`, `crz`, `cp` to $3 \times 10^{-7}$, the larger figure reflecting error accumulated over the two or three RZ rotations each contains.

Note that `ct` and `ctdg` are *exact* gates that nonetheless take the approximate path, because controlled-$T$ requires an arbitrary rotation. If a $\pi/8$-controlled phase is needed cheaply, express it directly rather than through `ct`.

### Verification

Every gate in both tables above is checked by `scripts/260829_verify_gate_decompositions.py`, which transpiles each one and compares the resulting unitary against Qiskit's reference for the same gate, up to global phase. Exact lowerings are required to agree to $10^{-9}$ and all land at machine precision; synthesized ones are required to agree to $10^{-5}$ and land between $7 \times 10^{-11}$ for a single rotation and $3 \times 10^{-7}$ for gates such as `u3` and `cu` that accumulate error over two or three.

This check exists because two defects went unnoticed without it, both fixed on 2026-08-29 and both now covered by explicit regression guards:

- **`cs` lowered to the wrong gate.** It emitted `s; cx; sdg; cx; s`, which is $\mathrm{diag}(1,1,1,-1) = \mathrm{CZ}$ rather than $\mathrm{diag}(1,1,1,i) = \mathrm{CS}$: the phases were $\pi/2$ where they should have been $\pi/4$. The failure was silent, producing a valid Clifford+T circuit computing the wrong unitary at zero $T$ cost. Now `t; cx; tdg; cx; t`, mirroring `csdg`, at 3 $T$.
- **`ecr` did not lower at all.** It had no case in the decomposition pass and passed through unchanged, so `to_clifford_t` returned a circuit still containing `ecr` while reporting success. Now lowered exactly.

Not Supported from `qelib1.inc`
-------------------------------

Five names in the standard library have no builtin here, and a file using them fails with `Unknown gate`:

| name | note |
|---|---|
| `u0` | idle/identity with a length parameter |
| `rc3x` | relative-phase 3-controlled X |
| `c3x`, `c3sqrtx`, `c4x` | fixed-size multi-controlled gates |

**Workaround.** Because the parser inlines user `gate` definitions, a file that carries the `qelib1.inc` definition of any of these parses normally. Qiskit's QASM 2 exporter emits such definitions in the preamble for non-standard gates, so exported files generally work without modification. Verified: a file defining `u0` and `rc3x` inline parses and inlines to 19 operations over `u1`, `u2`, `u3` and `cx`.

For multi-controlled X specifically, prefer `mcx` — it is variable-arity, and it is the form the lowering strategies understand.

Internal Representation Notes
-----------------------------

Two IR details are visible in output or in `count_ops()` but are not part of the input dialect.

- **`P4`, `P8`, `P16`** are internal $\pi/4$, $\pi/8$ and $\pi/16$ rotations produced by the Clifford reduction pass. They are **not parseable** — `p4 q[0];` fails with `Unknown gate` — but they are printed as ordinary OpenQASM 2: `t`, `tdg`, or `rx(pi/4)`, `rz(-pi/8)`, `rx(pi/16)` and so on, according to their axis and dagger flags. Emitted files therefore stay valid QASM 2. Round-tripping is semantically faithful but lossy in one respect: `rx(pi/4)` reads back as a generic `RX`, not as a `P4`, so the T-class bucketing in `count_ops()` is not reproduced. Note that `count_ops()` reports these under the key `rx(pi/4)` rather than `p4`.
- **`SWAP_BASIS`** is declared in the `Operation::Type` enum but is not produced or consumed anywhere in the codebase. It is reserved and currently unused; do not rely on it.

Interoperability Summary
------------------------

| you want | do this |
|---|---|
| a file any OpenQASM 2 tool can read | avoid `mcx`, `cs`/`csdg`/`ct`/`ctdg`, `ecr`, `ryy`, and any conditional |
| feed-forward / dynamic circuits | accept the OpenQASM 3 `if` syntax; the marker comment identifies such files |
| PBC output | it is NWQEC-only and does not round-trip through other tools |
| to import a Qiskit-exported QASM 2 file | generally works as-is, since Qiskit emits `gate` definitions for anything outside the builtin set |
| multi-controlled X | use `mcx`, not `c3x`/`c4x` |
