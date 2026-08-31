#pragma once

#include "pass_template.hpp"
#include "nwqec/tableau/vtab.hpp"
#include "nwqec/core/pauli_op.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <optional>
#include <unordered_map>

namespace NWQEC
{
    class PbcPass : public Pass
    {
    private:
        bool keep_cx;

    public:
        PbcPass(bool keep_cx = false) : keep_cx(keep_cx) {}

        /**
         * @brief Per-row metadata for an emitted PBC operation.
         *
         * The tableau only stores Pauli rows; everything else about what that row
         * becomes on output lives here, in the same insertion order.
         */
        struct StabMeta
        {
            Operation::Type type = Operation::Type::T_PAULI;
            std::vector<size_t> bits;                       // measurement destination
            std::optional<ClassicalCondition> condition;    // guard, if any
        };

        bool run(Circuit &circuit) override
        {
            std::vector<Operation> prepared = prepare_operations(circuit);

            std::vector<const Operation *> pbc_operations;
            pbc_operations.reserve(prepared.size());
            for (const auto &op : prepared)
            {
                if (op.get_type() == Operation::Type::T_PAULI ||
                    op.get_type() == Operation::Type::S_PAULI ||
                    op.get_type() == Operation::Type::M_PAULI ||
                    op.get_type() == Operation::Type::Z_PAULI)
                {
                    return false; // Already a PBC circuit
                }

                // Cancelling a self-inverse pair is only valid when neither is guarded
                // and nothing guarded sits between them.
                if (!pbc_operations.empty() &&
                    !op.is_conditional() && !pbc_operations.back()->is_conditional() &&
                    is_self_inverse_gate(op.get_type()) &&
                    same_gate_target(*pbc_operations.back(), op))
                {
                    pbc_operations.pop_back();
                    continue;
                }
                pbc_operations.push_back(&op);
            }

            size_t n_qubits = circuit.get_num_qubits();
            size_t n_gate_stabs = 0;

            for (const Operation *op_ptr : pbc_operations)
            {
                const Operation &op = *op_ptr;
                if (op.get_type() == Operation::Type::CCX)
                    n_gate_stabs += 7;
                else if (keep_cx && op.get_type() == Operation::Type::CX)
                    ++n_gate_stabs;
                else if (op.get_type() == Operation::Type::T || op.get_type() == Operation::Type::TDG)
                    ++n_gate_stabs;
                else if (op.get_type() == Operation::Type::MEASURE)
                    n_gate_stabs += op.get_qubits().size();
                else if (op.is_conditional())
                    n_gate_stabs += conditional_row_count(op);
            }

            VTab tableau(n_qubits, n_gate_stabs);
            std::vector<StabMeta> meta;
            meta.reserve(n_gate_stabs);

            for (auto it = pbc_operations.rbegin(); it != pbc_operations.rend(); ++it)
            {
                const Operation &op = **it;

                if (op.get_type() == Operation::Type::BARRIER)
                {
                    continue;
                }

                // A guarded operation cannot join the Clifford frame: whether it acts
                // depends on a measurement outcome unknown at compile time. It is
                // emitted in place instead, with its Pauli expressed in the current
                // frame, and deliberately NOT applied to the tableau.
                if (op.is_conditional())
                {
                    add_conditional_rows(op, n_qubits, tableau, meta);
                    continue;
                }

                if (op.get_type() == Operation::Type::MEASURE)
                {
                    // Mid-circuit measurement becomes an in-place Pauli measurement,
                    // its basis carried by the frame like any other row.
                    const auto &qubits = op.get_qubits();
                    const auto &bits = op.get_bits();
                    for (size_t i = qubits.size(); i-- > 0;)
                    {
                        PauliOp stab(n_qubits);
                        stab.set_r(false);
                        stab.add_z(qubits[i]);
                        tableau.add_stab(stab);
                        StabMeta m;
                        m.type = Operation::Type::M_PAULI;
                        if (i < bits.size())
                            m.bits = {bits[i]};
                        meta.push_back(m);
                    }
                    continue;
                }

                if (op.get_type() == Operation::Type::RESET)
                {
                    throw std::runtime_error(
                        "PbcPass: RESET on qubit " + std::to_string(op.get_qubits()[0]) +
                        " could not be rewritten as a measurement-conditioned X. RESET is not "
                        "unitary, so no Clifford frame can be swept through it. Ensure each RESET "
                        "directly follows a MEASURE of the same qubit.");
                }

                if (op.get_type() == Operation::Type::CCX)
                {
                    const auto &qubits = op.get_qubits();
                    tableau.add_ccx_stabs(qubits[0], qubits[1], qubits[2]);
                    for (size_t i = 0; i < 7; ++i)
                        meta.push_back(StabMeta{Operation::Type::T_PAULI, {}, std::nullopt});
                }
                else if (keep_cx && op.get_type() == Operation::Type::CX)
                {
                    const auto &qubits = op.get_qubits();
                    tableau.apply_clifford_gate(Operation::Type::SDG, qubits[0]);
                    tableau.apply_clifford_gate(Operation::Type::SXDG, qubits[1]);

                    PauliOp stab(n_qubits);
                    stab.set_r(false); // S_PAULI stabilizer
                    stab.add_z(qubits[0]);
                    stab.add_x(qubits[1]);
                    tableau.add_stab(stab);
                    meta.push_back(StabMeta{Operation::Type::S_PAULI, {}, std::nullopt});
                }
                else
                {
                    const auto &qubits = op.get_qubits();
                    if (op.get_type() == Operation::Type::T || op.get_type() == Operation::Type::TDG)
                    {
                        uint8_t phase = (op.get_type() == Operation::Type::T) ? 0 : 1;
                        tableau.add_t_stab(qubits[0], phase);
                        meta.push_back(StabMeta{Operation::Type::T_PAULI, {}, std::nullopt});
                    }
                    else
                    {
                        tableau.apply_clifford_gate(op.get_type(), qubits[0], qubits.size() > 1 ? qubits[1] : SIZE_MAX);
                    }
                }
            }
            std::vector<PauliOp> stabilizers = tableau.get_paili_ops();

            update_circuit(stabilizers, circuit, meta);
            return true;
        }

        std::string get_name() const override
        {
            return "PBC Pass";
        }

    private:
        /**
         * @brief Rewrite RESET into a measurement-conditioned X, and copy the rest.
         *
         * RESET is not unitary, so no Clifford frame can be swept through it. When it
         * directly follows a measurement of the same qubit -- which is what the Jones
         * gadget always produces -- it is exactly `if (m == 1) X(a)`, which is a
         * conditional Clifford and is handled like any other guarded operation.
         */
        std::vector<Operation> prepare_operations(const Circuit &circuit) const
        {
            const auto &ops = circuit.get_operations();
            std::vector<Operation> out;
            out.reserve(ops.size());

            // Last measurement destination bit seen for each qubit, if that qubit has
            // not been touched since.
            std::unordered_map<size_t, size_t> pending_measure;

            for (const auto &op : ops)
            {
                if (op.get_type() == Operation::Type::RESET && !op.is_conditional())
                {
                    size_t q = op.get_qubits()[0];
                    auto it = pending_measure.find(q);
                    if (it != pending_measure.end())
                    {
                        Operation x(Operation::Type::X, {q});
                        x.set_condition(ClassicalCondition({it->second}, 1));
                        out.push_back(x);
                        pending_measure.erase(it);
                        continue;
                    }
                    out.push_back(op); // left for run() to reject with a clear message
                    continue;
                }

                if (op.get_type() == Operation::Type::MEASURE)
                {
                    const auto &qs = op.get_qubits();
                    const auto &bs = op.get_bits();
                    for (size_t i = 0; i < qs.size() && i < bs.size(); ++i)
                        pending_measure[qs[i]] = bs[i];
                }
                else
                {
                    for (size_t q : op.get_qubits())
                        pending_measure.erase(q);
                }
                out.push_back(op);
            }
            return out;
        }

        /**
         * @brief How many Pauli rows a guarded Clifford expands into.
         */
        size_t conditional_row_count(const Operation &op) const
        {
            switch (op.get_type())
            {
            case Operation::Type::CZ:
            case Operation::Type::CX:
            case Operation::Type::H:
                return 3;
            case Operation::Type::X:
            case Operation::Type::Y:
            case Operation::Type::Z:
            case Operation::Type::S:
            case Operation::Type::SDG:
                return 1;
            default:
                throw std::runtime_error(
                    "PbcPass: no Pauli-rotation form is defined for a classically conditioned '" +
                    op.get_type_name() + "'. Supported: H, CX, CZ, X, Y, Z, S, SDG.");
            }
        }

        /**
         * @brief Emit a guarded Clifford as guarded Pauli rotations.
         *
         * Writing R(P,t) for a rotation by angle t about Pauli string P, and reading
         * products right-to-left as circuits:
         *
         *   CZ(a,b) = R(Z_a, pi/2) R(Z_b, pi/2) R(Z_a Z_b, -pi/2)
         *   CX(a,b) = R(Z_a, pi/2) R(X_b, pi/2) R(Z_a X_b, -pi/2)
         *   H(a)    = R(Z_a, pi/2) R(X_a, pi/2) R(Z_a, pi/2)
         *   X(a)    = R(X_a, pi)
         *
         * each up to a global phase, which is irrelevant to the measurement statistics
         * this pass is verified against. A positive angle is r = false and a negative
         * angle r = true, matching how T and TDG are encoded.
         *
         * Rows are pushed here in reverse circuit order, because the caller walks the
         * circuit backwards and update_circuit reverses once more; that is why the
         * leftmost factor above is pushed first. Order matters for CX and H, whose
         * factors do not commute.
         *
         * Every emitted row carries the original guard.
         */
        void add_conditional_rows(const Operation &op, size_t n_qubits,
                                  VTab &tableau, std::vector<StabMeta> &meta) const
        {
            const auto &q = op.get_qubits();
            const auto &cond = *op.get_condition();

            auto push = [&](Operation::Type type, const PauliOp &p)
            {
                tableau.add_stab(p);
                meta.push_back(StabMeta{type, {}, cond});
            };

            switch (op.get_type())
            {
            case Operation::Type::CZ:
            {
                PauliOp za(n_qubits); za.set_r(false); za.add_z(q[0]);
                PauliOp zb(n_qubits); zb.set_r(false); zb.add_z(q[1]);
                PauliOp zz(n_qubits); zz.set_r(true);  zz.add_z(q[0]); zz.add_z(q[1]);
                push(Operation::Type::S_PAULI, za);
                push(Operation::Type::S_PAULI, zb);
                push(Operation::Type::S_PAULI, zz);
                break;
            }
            case Operation::Type::CX:
            {
                PauliOp za(n_qubits); za.set_r(false); za.add_z(q[0]);
                PauliOp xb(n_qubits); xb.set_r(false); xb.add_x(q[1]);
                PauliOp zx(n_qubits); zx.set_r(true);  zx.add_z(q[0]); zx.add_x(q[1]);
                push(Operation::Type::S_PAULI, za);
                push(Operation::Type::S_PAULI, xb);
                push(Operation::Type::S_PAULI, zx);
                break;
            }
            case Operation::Type::H:
            {
                PauliOp z1(n_qubits); z1.set_r(false); z1.add_z(q[0]);
                PauliOp xm(n_qubits); xm.set_r(false); xm.add_x(q[0]);
                PauliOp z2(n_qubits); z2.set_r(false); z2.add_z(q[0]);
                push(Operation::Type::S_PAULI, z1);
                push(Operation::Type::S_PAULI, xm);
                push(Operation::Type::S_PAULI, z2);
                break;
            }
            case Operation::Type::X:
            {
                PauliOp p(n_qubits); p.set_r(false); p.add_x(q[0]);
                push(Operation::Type::Z_PAULI, p);
                break;
            }
            case Operation::Type::Y:
            {
                PauliOp p(n_qubits); p.set_r(false); p.add_x(q[0]); p.add_z(q[0]);
                push(Operation::Type::Z_PAULI, p);
                break;
            }
            case Operation::Type::Z:
            {
                PauliOp p(n_qubits); p.set_r(false); p.add_z(q[0]);
                push(Operation::Type::Z_PAULI, p);
                break;
            }
            case Operation::Type::S:
            case Operation::Type::SDG:
            {
                PauliOp p(n_qubits);
                p.set_r(op.get_type() == Operation::Type::SDG);
                p.add_z(q[0]);
                push(Operation::Type::S_PAULI, p);
                break;
            }
            default:
                throw std::runtime_error(
                    "PbcPass: no Pauli-rotation form is defined for a classically conditioned '" +
                    op.get_type_name() + "'.");
            }
        }

        bool is_self_inverse_gate(Operation::Type type) const
        {
            return type == Operation::Type::H ||
                   type == Operation::Type::X ||
                   type == Operation::Type::Y ||
                   type == Operation::Type::Z ||
                   type == Operation::Type::CX;
        }

        bool same_gate_target(const Operation &lhs, const Operation &rhs) const
        {
            return lhs.get_type() == rhs.get_type() &&
                   lhs.get_qubits() == rhs.get_qubits();
        }

        void update_circuit(std::vector<PauliOp> &stabilizers, Circuit &circuit,
                            std::vector<StabMeta> &meta)
        {
            size_t n_qubits = circuit.get_num_qubits();
            size_t n_bits = circuit.get_num_bits();

            Circuit new_circuit;
            new_circuit.add_qreg("q", n_qubits);
            if (n_bits > 0)
                new_circuit.add_creg("c", n_bits);

            // Rows were inserted walking the circuit backwards, so reversing here
            // restores circuit order.
            std::reverse(meta.begin(), meta.end());
            size_t stab_idx = 0;

            assert(stabilizers.size() == n_qubits + meta.size());
            for (size_t pos = stabilizers.size(); pos-- > n_qubits;)
            {
                const PauliOp &pauli_op = stabilizers[pos];
                const StabMeta &m = meta[stab_idx++];
                Operation op(m.type, {}, {}, m.bits, pauli_op);
                if (m.condition.has_value())
                    op.set_condition(*m.condition);
                new_circuit.add_operation(op);
            }

            // Terminal Pauli measurements, which carry the leftover Clifford frame.
            for (size_t i = 0; i < n_qubits; ++i)
            {
                new_circuit.add_operation(Operation(Operation::Type::M_PAULI, {}, {}, {}, stabilizers[i]));
            }

            circuit = std::move(new_circuit);
        }
    };

} // namespace NWQEC
