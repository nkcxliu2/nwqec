#pragma once

#include "nwqec/core/operation.hpp"
#include "nwqec/core/classical_condition.hpp"

#include <vector>

namespace NWQEC
{
    /**
     * @brief Gidney's temporary logical AND.
     *
     * Compute writes w = a AND b into a clean ancilla w for 4 T; erasure measures the
     * ancilla in the X basis and repairs the leftover phase with a classically
     * conditioned CZ, for 0 T. Erasure is only valid when w is a genuine clean ancilla
     * holding nothing but this AND: it destroys w's superposition, so it must never be
     * applied to a live data qubit.
     *
     * Craig Gidney, "Halving the cost of quantum addition", Quantum 2, 74 (2018).
     */
    inline std::vector<Operation> and_compute(size_t a, size_t b, size_t w)
    {
        return {Operation(Operation::Type::RCCX, {a, b, w}),
                Operation(Operation::Type::SDG, {w})};
    }

    inline std::vector<Operation> and_uncompute(size_t a, size_t b, size_t w,
                                                size_t result_bit, bool emit_reset = true)
    {
        std::vector<Operation> ops;
        ops.push_back(Operation(Operation::Type::H, {w}));
        ops.push_back(Operation(Operation::Type::MEASURE, {w}, {}, {result_bit}));
        Operation fix(Operation::Type::CZ, {a, b});
        fix.set_condition(ClassicalCondition({result_bit}, 1));
        ops.push_back(std::move(fix));
        if (emit_reset)
        {
            ops.push_back(Operation(Operation::Type::RESET, {w}));
        }
        return ops;
    }

} // namespace NWQEC
