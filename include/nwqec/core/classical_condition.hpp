#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace NWQEC
{
    /**
     * @brief A condition on classical bits guarding the execution of an operation.
     *
     * `bits` holds global classical bit indices, least-significant first: bit i of
     * `value` is compared against the measured value of bits[i]. The operation
     * executes iff every comparison holds.
     *
     * A single-bit condition is {bits = {m}, value = 1}. A whole-register
     * comparison `if (c == k)` is {bits = [start .. start+size-1], value = k}.
     */
    struct ClassicalCondition
    {
        static constexpr size_t MAX_BITS = 64;

        std::vector<size_t> bits;
        uint64_t value = 0;

        ClassicalCondition() = default;

        ClassicalCondition(std::vector<size_t> condition_bits, uint64_t condition_value)
            : bits(std::move(condition_bits)), value(condition_value)
        {
            validate();
        }

        void validate() const
        {
            if (bits.empty())
            {
                throw std::invalid_argument("ClassicalCondition must reference at least one classical bit");
            }
            if (bits.size() > MAX_BITS)
            {
                throw std::invalid_argument("ClassicalCondition supports at most 64 classical bits");
            }
            if (bits.size() < MAX_BITS && value >= (uint64_t(1) << bits.size()))
            {
                throw std::invalid_argument("ClassicalCondition value does not fit in the referenced bit width");
            }
        }

        bool operator==(const ClassicalCondition &other) const
        {
            return bits == other.bits && value == other.value;
        }

        bool operator!=(const ClassicalCondition &other) const { return !(*this == other); }
    };

} // namespace NWQEC
