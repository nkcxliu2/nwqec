#pragma once

#include "operation.hpp"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <unordered_map>
#include <algorithm>

namespace NWQEC
{
    /**
     * Statistics for basis operations during depth calculation
     */
    struct BasisStatistics
    {
        size_t total_operations = 0;
        size_t z_basis_operations = 0;
        size_t x_basis_operations = 0;
        size_t y_operations = 0;
        size_t basis_changes = 0;

        void clear()
        {
            total_operations = 0;
            z_basis_operations = 0;
            x_basis_operations = 0;
            y_operations = 0;
            basis_changes = 0;
        }
    };

    /**
     * Represents a flattened quantum circuit with only elementary gates
     */
    class Circuit
    {
    private:
        size_t num_qubits = 0;
        size_t num_bits = 0;
        std::vector<Operation> operations;

        // Maps register names to their starting indices
        std::map<std::string, size_t> qubit_register_map;
        std::map<std::string, size_t> bit_register_map;

        // Maps register names to their sizes
        std::map<std::string, size_t> qubit_reg_size_map;
        std::map<std::string, size_t> bit_reg_size_map;

        // Maps user-defined gate names to their definitions
        std::map<std::string, std::vector<Operation>> gate_definitions;

        // Track if circuit contains only Clifford+T gates
        mutable bool is_clifford_t_circuit = true;

        // Track if circuit contains any classically conditioned operation
        bool has_feedforward_ops = false;

        // Basis statistics from last basis_aware_depth calculation
        mutable BasisStatistics basis_stats;

    public:
        std::vector<double> distinct_rz_angles;
        std::map<size_t, size_t> rz_angle_map;

        Circuit() = default;
        virtual ~Circuit() = default; // Good practice to have a virtual function

        // Register management
        virtual void add_qreg(const std::string &name, size_t size)
        {
            qubit_register_map[name] = num_qubits;
            qubit_reg_size_map[name] = size;
            num_qubits += size;
        }

        virtual void add_creg(const std::string &name, size_t size)
        {
            bit_register_map[name] = num_bits;
            bit_reg_size_map[name] = size;
            num_bits += size;
        }

        // Get global qubit index from register name and local index
        size_t get_qubit_index(const std::string &reg_name, size_t local_index) const
        {
            auto it = qubit_register_map.find(reg_name);
            if (it == qubit_register_map.end())
            {
                throw std::runtime_error("Unknown quantum register: " + reg_name);
            }
            size_t startIdx = it->second;
            return startIdx + local_index;
        }

        size_t get_qubit_reg_size(const std::string &reg_name) const
        {
            auto it = qubit_reg_size_map.find(reg_name);
            if (it == qubit_reg_size_map.end())
            {
                throw std::runtime_error("Unknown quantum register: " + reg_name);
            }
            return it->second;
        }

        // Get global bit index from register name and local index
        size_t get_bit_index(const std::string &reg_name, size_t local_index) const
        {
            auto it = bit_register_map.find(reg_name);
            if (it == bit_register_map.end())
            {
                throw std::runtime_error("Unknown classical register: " + reg_name);
            }
            size_t startIdx = it->second;
            return startIdx + local_index;
        }

        size_t get_bit_reg_size(const std::string &reg_name) const
        {
            auto it = bit_reg_size_map.find(reg_name);
            if (it == bit_reg_size_map.end())
            {
                throw std::runtime_error("Unknown classical register: " + reg_name);
            }
            return it->second;
        }

        // Gate definition management
        void define_gate(const std::string &name, const std::vector<Operation> &definition)
        {
            gate_definitions[name] = definition;
        }

        bool is_user_defined_gate(const std::string &name) const
        {
            return gate_definitions.find(name) != gate_definitions.end();
        }

        // Expand a user-defined gate with actual qubits
        void expand_gate(const std::string &name, const std::vector<size_t> &actual_qubits)
        {
            auto it = gate_definitions.find(name);
            if (it == gate_definitions.end())
            {
                throw std::runtime_error("Unknown user-defined gate: " + name);
            }

            const std::vector<Operation> &definition = it->second;

            // Apply the gate definition with the specified qubits
            for (const auto &op : definition)
            {
                std::vector<size_t> mapped_qubits;
                for (size_t formalIdx : op.get_qubits())
                {
                    mapped_qubits.push_back(actual_qubits[formalIdx]);
                }

                add_operation(Operation(op.get_type(), mapped_qubits, op.get_parameters(), op.get_bits(),
                                        op.get_pauli_op(), op.get_dagger(), op.get_x_rotation(),
                                        op.get_condition()));
            }
        }

        // Add operations to the circuit
        void reserve_operations(size_t count)
        {
            operations.reserve(count);
        }

        virtual void add_operation(Operation operation)
        {
            // Update num_qubits and num_bits if new operation extends the circuit
            for (size_t q_idx : operation.get_qubits())
            {
                if (q_idx >= num_qubits)
                {
                    num_qubits = q_idx + 1;
                }
            }
            for (size_t b_idx : operation.get_bits())
            {
                if (b_idx >= num_bits)
                {
                    num_bits = b_idx + 1;
                }
            }
            // Condition bits are read, not written, but still occupy classical storage
            // and must be covered by the declared creg.
            if (operation.is_conditional())
            {
                for (size_t b_idx : operation.get_condition()->bits)
                {
                    if (b_idx >= num_bits)
                    {
                        num_bits = b_idx + 1;
                    }
                }
                has_feedforward_ops = true;
            }

            // Update Clifford+T tracking
            if (is_clifford_t_circuit && !is_clifford_t_operation(operation.get_type()))
            {
                is_clifford_t_circuit = false;
            }

            operations.push_back(std::move(operation));
        }

        /**
         * @brief Apply a classical condition to every operation from begin_idx onward.
         *
         * Used by the parser to guard the body of an `if` statement after it has been
         * parsed through the ordinary gate path, and by lowering passes that emit a
         * conditional gadget. Throws if any operation in the range cannot be guarded.
         */
        void set_condition_on_operations(size_t begin_idx, const ClassicalCondition &condition)
        {
            condition.validate();
            for (size_t i = begin_idx; i < operations.size(); ++i)
            {
                operations[i].set_condition(condition);
            }
            if (begin_idx < operations.size())
            {
                has_feedforward_ops = true;
                for (size_t b_idx : condition.bits)
                {
                    if (b_idx >= num_bits)
                    {
                        num_bits = b_idx + 1;
                    }
                }
            }
        }

        // Allow derived classes to set the full list of operations.
        // Useful for transformations that rebuild the operation list.
        void set_operations_list(std::vector<Operation> new_ops)
        {
            operations = std::move(new_ops);
            // Recalculate num_qubits and num_bits based on the new operations list
            num_qubits = 0;
            num_bits = 0;
            is_clifford_t_circuit = true;
            has_feedforward_ops = false;

            for (const auto &op : operations)
            {
                for (size_t q_idx : op.get_qubits())
                {
                    if (q_idx >= num_qubits)
                    {
                        num_qubits = q_idx + 1;
                    }
                }
                for (size_t b_idx : op.get_bits())
                {
                    if (b_idx >= num_bits)
                    {
                        num_bits = b_idx + 1;
                    }
                }
                if (op.is_conditional())
                {
                    for (size_t b_idx : op.get_condition()->bits)
                    {
                        if (b_idx >= num_bits)
                        {
                            num_bits = b_idx + 1;
                        }
                    }
                    has_feedforward_ops = true;
                }

                // Check Clifford+T status
                if (!is_clifford_t_operation(op.get_type()))
                {
                    is_clifford_t_circuit = false;
                }
            }
        }
        // Helper method to check if an operation is a Clifford+T gate
        bool is_clifford_t_operation(Operation::Type type) const
        {
            switch (type)
            {
            // Clifford gates
            case Operation::Type::H:
            case Operation::Type::S:
            case Operation::Type::SDG:
            case Operation::Type::X:
            case Operation::Type::Y:
            case Operation::Type::Z:
            case Operation::Type::SX:
            case Operation::Type::SXDG:
            case Operation::Type::CX:
            // T gates
            case Operation::Type::T:
            case Operation::Type::TDG:
            case Operation::Type::P4:
            case Operation::Type::T_PAULI:
            case Operation::Type::S_PAULI:
            // Ignored operations
            case Operation::Type::MEASURE:
            case Operation::Type::RESET:
            case Operation::Type::BARRIER:
            case Operation::Type::M_PAULI:
                return true;
            default:
                return false;
            }
        }

        // Getters
        size_t get_num_qubits() const { return num_qubits; }
        size_t get_num_bits() const { return num_bits; }
        const std::vector<Operation> &get_operations() const { return operations; }

        // Check if the circuit contains only Clifford+T gates
        bool is_clifford_t() const { return is_clifford_t_circuit; }

        // ---- dynamic-circuit queries ----
        // True if any operation is classically conditioned (measurement feed-forward).
        bool has_feedforward() const { return has_feedforward_ops; }

        bool has_measurement() const
        {
            for (const auto &op : operations)
            {
                if (op.is_measurement())
                    return true;
            }
            return false;
        }

        bool is_dynamic() const { return has_feedforward_ops || has_measurement(); }

        // Get count for a specific operation type
        size_t get_operation_count(Operation::Type type) const
        {
            size_t count = 0;
            for (const auto &op : operations)
            {
                if (op.get_type() == type)
                {
                    count++;
                }
            }
            return count;
        }

        // Update the qubit and bit counts after removing unused resources
        void update_qubit_and_bit_counts(size_t new_num_qubits, size_t new_num_bits)
        {
            num_qubits = new_num_qubits;
            num_bits = new_num_bits;
        }

        // ===================== Circuit Analysis Methods =====================
        // These methods are implemented in src/core/circuit_analysis.cpp

        /**
         * @brief Calculate circuit depth (critical path length)
         * @return Circuit depth in terms of gate layers
         */
        size_t depth() const
        {
            std::vector<size_t> depth_counts(num_qubits, 0); // track depth of each qubit

            for (const auto &operation : operations)
            {
                const auto &qubits = operation.get_qubits();
                size_t current_depth = 0;

                for (size_t qubit : qubits)
                {
                    if (qubit < depth_counts.size())
                        current_depth = std::max(depth_counts[qubit], current_depth);
                }
                for (size_t qubit : qubits)
                {
                    if (qubit < depth_counts.size())
                        depth_counts[qubit] = current_depth + 1;
                }
            }

            size_t max_depth = 0;
            for (size_t depth : depth_counts)
            {
                max_depth = std::max(max_depth, depth);
            }

            return max_depth;
        }

        // Helper function to get gate duration for different operation types
        double get_gate_duration(Operation::Type type, double code_distance) const
        {
            switch (type)
            {
            case Operation::Type::CX:
                return 3 * code_distance + 4;
            case Operation::Type::H:
                return 3 * code_distance + 4;
            case Operation::Type::S:
                return 1.5 * code_distance + 3;
            case Operation::Type::SDG:
                return 1.5 * code_distance + 3;
            case Operation::Type::SX:
                return 1.5 * code_distance + 3;
            case Operation::Type::SXDG:
                return 1.5 * code_distance + 3;
            case Operation::Type::T:
                return 2.5 * code_distance + 4;
            case Operation::Type::TDG:
                return 2.5 * code_distance + 4;
            case Operation::Type::P4:
                return 2.5 * code_distance + 4;
            default:
                return code_distance; // Default duration for all other gates
            }
        }

        /**
         * @brief Calculate circuit execution duration based on gate timings
         * @return Total execution duration in time units
         */
        double duration(double code_distance) const
        {
            std::vector<double> duration_counts(num_qubits, 0.0); // track duration of each qubit

            for (const auto &operation : operations)
            {
                const auto &qubits = operation.get_qubits();

                // Get gate duration based on type
                double gate_duration = get_gate_duration(operation.get_type(), code_distance);

                double current_duration = 0.0;

                for (size_t qubit : qubits)
                {
                    if (qubit < duration_counts.size())
                        current_duration = std::max(duration_counts[qubit], current_duration);
                }
                for (size_t qubit : qubits)
                {
                    if (qubit < duration_counts.size())
                        duration_counts[qubit] = current_duration + gate_duration;
                }
            }

            double max_duration = 0.0;
            for (double duration : duration_counts)
            {
                max_duration = std::max(max_duration, duration);
            }

            return max_duration;
        }

        /**
         * @brief Print comprehensive circuit statistics
         * @param os Output stream to print to
         * @param skip_pauli Whether to skip Pauli gates in calculations
         */
        void print_stats(std::ostream &os = std::cout) const
        {
            // Get operation counts using the count_ops method
            auto op_counts = count_ops();

            os << "==================================================\n";
            os << "Circuit Statistics\n";
            os << "==================================================\n";

            os << "Basic Circuit Information:\n";
            os << "  Number of qubits: " << num_qubits << "\n";
            os << "  Number of classical bits: " << num_bits << "\n";
            os << "  Total gates: " << operations.size() << "\n";
            os << "  Circuit depth: " << depth() << " gates\n";
            os << "  Clifford+T circuit: " << (is_clifford_t_circuit ? "Yes" : "No") << "\n";

            os << "\nGate Count Breakdown:\n";

            // Priority order for gate operations
            const std::vector<std::string> priority_ops = {
                "t", "tdg", "rx(pi/4)", "rx(-pi/4)", "t_pauli", "s", "sdg", "s_pauli", "h", "cx", "ccx", "measure", "m_pauli"};

            // First print priority operations in order
            for (const auto &op_name : priority_ops)
            {
                if (op_counts.find(op_name) != op_counts.end() && op_counts[op_name] > 0)
                {
                    os << "    " << op_name << ": " << op_counts[op_name] << "\n";
                    op_counts.erase(op_name); // Remove printed op from the map
                }
            }

            // Then print remaining operations
            for (const auto &[op_name, count] : op_counts)
            {
                os << "    " << op_name << ": " << count << "\n";
            }

            // Print derived statistics
            size_t t_gates = get_operation_count(Operation::Type::T) +
                             get_operation_count(Operation::Type::TDG) +
                             get_operation_count(Operation::Type::T_PAULI) +
                             get_operation_count(Operation::Type::P4);

            size_t two_qubit_gates = get_operation_count(Operation::Type::CX) +
                                     get_operation_count(Operation::Type::CY) +
                                     get_operation_count(Operation::Type::CZ) +
                                     get_operation_count(Operation::Type::SWAP) +
                                     get_operation_count(Operation::Type::ECR);

            os << "\nDerived Statistics:\n";
            os << "  Total T-type gates: " << t_gates << "\n";
            os << "  Total two-qubit gates: " << two_qubit_gates << "\n";

            if (operations.size() > 0)
            {
                double t_gate_ratio = static_cast<double>(t_gates) / operations.size() * 100;
                double two_qubit_ratio = static_cast<double>(two_qubit_gates) / operations.size() * 100;

                os << "  T-gate ratio: " << std::fixed << std::setprecision(1) << t_gate_ratio << "%\n";
                os << "  Two-qubit gate ratio: " << std::fixed << std::setprecision(1) << two_qubit_ratio << "%\n";
            }

            os << "==================================================\n\n";
        }

        /**
         * @brief Print circuit statistics with Clifford+T focus
         * @param os Output stream to print to
         */
        void print_stats_ct(std::ostream &os = std::cout) const
        {
            // Get operation counts using the count_ops method
            auto op_counts = count_ops();

            os << "==================================================\n";
            os << "Clifford+T Circuit Statistics\n";
            os << "==================================================\n";

            os << "  Number of qubits: " << num_qubits << "\n";
            os << "  Total gates: " << operations.size() << "\n";
            os << "  Circuit depth: " << depth() << " gates\n";

            // Define gate categories with their gates
            struct GateCategory
            {
                std::string name;
                std::vector<std::string> gates;
            };

            std::vector<GateCategory> categories = {
                {"H", {"h"}},
                {"S", {"s", "sdg", "sx", "sxdg"}},
                {"T", {"t", "tdg", "rx(pi/4)", "rx(-pi/4)", "t_pauli"}},
                {"CX", {"cx"}},
                {"CCX", {"ccx"}},
                {"Pauli", {"x", "y", "z"}}};

            os << "\nGate Categories:\n";
            // Print each category
            for (const auto &category : categories)
            {
                size_t total = 0;

                // Calculate total for this category
                for (const auto &gate_name : category.gates)
                {
                    auto it = op_counts.find(gate_name);
                    if (it != op_counts.end() && it->second > 0)
                    {
                        total += it->second;
                        op_counts.erase(gate_name); // Remove printed gate from the map
                    }
                }

                os << "  " << std::setw(5) << category.name << ": " << total << "\n";
            }

            // Print remaining gates
            if (!op_counts.empty())
            {
                os << "\nOther gates:\n";
                for (const auto &[op_name, count] : op_counts)
                {
                    os << "    " << std::setw(10) << op_name << ": " << count << "\n";
                }
            }

            os << "==================================================\n\n";
        }

        /**
         * @brief Count operations by type and return as string-based map
         * @return Map of operation type names to their counts
         */
        std::unordered_map<std::string, size_t> count_ops() const
        {
            std::unordered_map<std::string, size_t> op_counts;

            // Count operations directly from the operations vector
            for (const auto &op : operations)
            {
                std::string op_name = op.get_type_name();
                op_counts[op_name]++;
            }

            return op_counts;
        }

        /**
         * @brief Render a classical condition as an OpenQASM 3 `if` header.
         *
         * A single-bit condition prints as `if (c[3] == 1)`. A condition covering the
         * whole flattened register prints as `if (c == 5)`. The indexed form is always
         * correct against the flattened register, which is why it is the default; the
         * whole-register form is only used when the condition genuinely spans every bit.
         */
        std::string condition_to_string(const ClassicalCondition &cond) const
        {
            if (cond.bits.size() == 1)
            {
                return "if (c[" + std::to_string(cond.bits[0]) + "] == " +
                       std::to_string(cond.value) + ")";
            }

            bool spans_register = (cond.bits.size() == num_bits);
            if (spans_register)
            {
                for (size_t i = 0; i < cond.bits.size(); ++i)
                {
                    if (cond.bits[i] != i)
                    {
                        spans_register = false;
                        break;
                    }
                }
            }
            if (spans_register)
            {
                return "if (c == " + std::to_string(cond.value) + ")";
            }

            // General case: conjunction over individually indexed bits.
            std::string expr = "if (";
            for (size_t i = 0; i < cond.bits.size(); ++i)
            {
                if (i > 0)
                    expr += " && ";
                expr += "c[" + std::to_string(cond.bits[i]) + "] == " +
                        std::to_string((cond.value >> i) & 1ULL);
            }
            return expr + ")";
        }

        // Print the circuit
        void print(std::ostream &os) const
        {
            // Print header. The header and include line are identical in both modes;
            // only conditional statements deviate from OpenQASM 2 (see plan section 10.3).
            os << "OPENQASM 2.0;\n";
            os << "include \"qelib1.inc\";\n";
            if (has_feedforward_ops)
            {
                os << "// NWQEC dynamic dialect: conditional statements use OpenQASM 3 "
                      "if-syntax.\n";
            }
            os << "\n";

            // Print registers
            if (num_qubits > 0)
            {
                os << "qreg q[" << num_qubits << "];\n";
            }
            if (num_bits > 0)
            {
                os << "creg c[" << num_bits << "];\n";
            }
            os << "\n";

            // Print operations. Adjacent operations sharing the identical condition are
            // emitted as a single block, so one logical conditional gate stays one
            // statement even after it is decomposed into many basis gates.
            for (size_t i = 0; i < operations.size();)
            {
                if (!operations[i].is_conditional())
                {
                    operations[i].print(os);
                    os << "\n";
                    ++i;
                    continue;
                }

                const ClassicalCondition &cond = *operations[i].get_condition();
                size_t j = i + 1;
                while (j < operations.size() && operations[j].is_conditional() &&
                       *operations[j].get_condition() == cond)
                {
                    ++j;
                }

                os << condition_to_string(cond) << " { ";
                for (size_t k = i; k < j; ++k)
                {
                    if (k > i)
                        os << " ";
                    operations[k].print(os);
                }
                os << " }\n";
                i = j;
            }
        }
    };

    /**
     * @brief Reject a circuit carrying classical feed-forward.
     *
     * Used by passes that operate only on unitary regions. Failing loudly is
     * deliberate: silently skipping a conditional operation would change the
     * circuit's semantics without any diagnostic.
     */
    inline void require_no_feedforward(const Circuit &circuit, const char *pass_name)
    {
        if (circuit.has_feedforward())
        {
            throw std::runtime_error(
                std::string(pass_name) +
                " does not support classically conditioned operations; "
                "this pipeline operates on unitary regions only.");
        }
    }

} // namespace NWQEC
