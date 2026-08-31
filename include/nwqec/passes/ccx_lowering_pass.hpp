#pragma once

#include "pass_template.hpp"
#include "nwqec/core/classical_condition.hpp"
#include "and_gadget.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace NWQEC
{
    /**
     * @brief How a high-level CCX is realized.
     */
    enum class CCXStrategy
    {
        Preserve,   // leave CCX untouched
        Standard7T, // current exact unitary lowering; 7 T, 6 CX, no ancilla
        Jones4T,    // one AND per Toffoli; 4 T each
        SharedAND   // one AND per group of Toffolis sharing a control pair; 4 T per group
    };

    struct CCXLoweringOptions
    {
        CCXStrategy strategy = CCXStrategy::Standard7T;
        // Whether an ancilla may be used at all. This is a property of the target, not
        // an optimization choice, so it overrides the strategy: with no ancilla the only
        // exact option is the ancilla-free seven-T decomposition.
        bool allow_ancilla = true;
        // Ancilla budget; 0 means unlimited, i.e. as many as the grouping requires.
        size_t max_ancillas = 0;
        bool emit_reset = true;
    };

    /**
     * @brief Lowers CCX by computing logical ANDs and reading them.
     *
     * A group is a maximal run of CCX operations that share a control pair {a,b} with
     * no intervening operation disturbing either control. For a group with ancilla w:
     *
     *     RCCX(a,b,w); SDG(w)                      compute w = a AND b     [4 T]
     *     CX(w, t_i)                               one per member          [0 T]
     *     H(w); MEASURE(w->m); if(m) CZ(a,b)       erase w                 [0 T]
     *     RESET(w)
     *
     * The substitution is exact because w holds a AND b, so CX(w,t) flips t exactly when
     * CCX(a,b,t) would. The cost is 4 T per group rather than per Toffoli.
     *
     * A group of size one is character for character the Jones gadget, so Jones4T is
     * simply this pass with every group split to size one.
     *
     * Control stability does double duty: it keeps a AND b valid at every later use, and
     * it makes the erasure's CZ(a,b) fixup correct, since the fixup is evaluated against
     * the same a and b that produced the AND.
     *
     * Craig Gidney, "Halving the cost of quantum addition", Quantum 2, 74 (2018).
     * Cody Jones, Phys. Rev. A 87, 022328 (2013), for the single-use case.
     */
    class CCXLoweringPass : public Pass
    {
    public:
        struct Group
        {
            size_t c0 = 0, c1 = 0;        // control pair, as first encountered
            std::vector<size_t> members;  // indices into the source operation list
            size_t ancilla = 0;
            size_t result_bit = 0;
            size_t begin() const { return members.front(); }
            size_t end() const { return members.back(); }
        };

        explicit CCXLoweringPass(CCXLoweringOptions options = {}) : options_(options) {}

        std::string get_name() const override { return "CCX Lowering Pass"; }

        size_t ancillas_added() const { return ancillas_added_; }
        size_t result_bits_added() const { return result_bits_added_; }
        const std::vector<Group> &groups() const { return groups_; }

        /**
         * @brief Qubits whose computational-basis value an operation changes.
         *
         * Diagonal gates and control roles are absent: they leave the value intact, so
         * they cannot break an AND. A classically conditioned operation is tested by the
         * gate it guards, but a guarded gate that could disturb a control is treated as
         * disturbing it whether or not it fires, since the outcome is unknown here.
         */
        static std::vector<size_t> disturbed_qubits(const Operation &op)
        {
            const auto &q = op.get_qubits();
            switch (op.get_type())
            {
            // Diagonal: the basis value is untouched.
            case Operation::Type::Z:
            case Operation::Type::S:
            case Operation::Type::SDG:
            case Operation::Type::T:
            case Operation::Type::TDG:
            case Operation::Type::RZ:
            case Operation::Type::P:
            case Operation::Type::U1:
            case Operation::Type::CZ:
            case Operation::Type::CS:
            case Operation::Type::CSDG:
            case Operation::Type::CT:
            case Operation::Type::CTDG:
            case Operation::Type::CP:
            case Operation::Type::CU1:
            case Operation::Type::CRZ:
            case Operation::Type::RZZ:
            case Operation::Type::ID:
                return {};
            // Controlled: only the target changes.
            case Operation::Type::CX:
            case Operation::Type::CY:
            case Operation::Type::CH:
            case Operation::Type::CSX:
            case Operation::Type::CRX:
            case Operation::Type::CRY:
            case Operation::Type::CU:
            case Operation::Type::CU3:
                return q.size() > 1 ? std::vector<size_t>{q[1]} : std::vector<size_t>{};
            case Operation::Type::CCX:
            case Operation::Type::RCCX:
                return q.size() > 2 ? std::vector<size_t>{q[2]} : std::vector<size_t>{};
            case Operation::Type::CSWAP:
                return q.size() > 2 ? std::vector<size_t>{q[1], q[2]} : std::vector<size_t>{};
            // Everything else, including X/Y/H/SX/RX/RY, SWAP, MEASURE, RESET, BARRIER
            // and the Pauli-frame operations, is treated as disturbing all its operands.
            default:
                return q;
            }
        }

        bool run(Circuit &circuit) override
        {
            if (options_.strategy != CCXStrategy::Jones4T &&
                options_.strategy != CCXStrategy::SharedAND)
            {
                return false; // Preserve and Standard7T are handled by DecomposePass
            }
            if (!options_.allow_ancilla)
            {
                return false; // no ancilla, so DecomposePass emits the seven-T circuit
            }

            const auto &ops = circuit.get_operations();
            groups_ = build_groups(ops);
            if (groups_.empty())
            {
                return false;
            }

            // Jones is this pass with every group split to size one.
            size_t budget = (options_.strategy == CCXStrategy::Jones4T)
                                ? 1
                                : options_.max_ancillas;
            if (budget == 1)
            {
                split_to_singletons(groups_);
            }
            else if (budget > 1)
            {
                enforce_budget(groups_, budget);
            }

            const size_t n_ancillas = assign_ancillas(groups_);
            const size_t nq_in = circuit.get_num_qubits();
            const size_t nb_in = circuit.get_num_bits();

            // Size the registers once: DAGCircuit rejects a qubit index beyond the
            // declared count, and every downstream pass rebuilds from these counts.
            Circuit out;
            out.add_qreg("q", nq_in + n_ancillas);
            out.add_creg("c", nb_in + n_ancillas);
            for (auto &g : groups_)
            {
                g.ancilla += nq_in;
                g.result_bit = nb_in + (g.ancilla - nq_in);
            }

            // member index -> (group, position within the group)
            std::vector<const Group *> first_of(ops.size(), nullptr);
            std::vector<const Group *> member_of(ops.size(), nullptr);
            std::vector<const Group *> last_of(ops.size(), nullptr);
            for (const auto &g : groups_)
            {
                for (size_t m : g.members)
                    member_of[m] = &g;
                first_of[g.begin()] = &g;
                last_of[g.end()] = &g;
            }

            for (size_t i = 0; i < ops.size(); ++i)
            {
                if (member_of[i] == nullptr)
                {
                    out.add_operation(ops[i]);
                    continue;
                }
                const Group &g = *member_of[i];
                if (first_of[i] == &g)
                {
                    emit_and_compute(g, out);
                }
                out.add_operation(Operation(Operation::Type::CX, {g.ancilla, ops[i].get_qubits()[2]}));
                if (last_of[i] == &g)
                {
                    emit_and_uncompute(g, out);
                }
            }

            ancillas_added_ = n_ancillas;
            result_bits_added_ = n_ancillas;
            circuit = std::move(out);
            return true;
        }

    private:
        /** Maximal runs of CCX sharing a control pair with both controls undisturbed. */
        std::vector<Group> build_groups(const std::vector<Operation> &ops) const
        {
            std::vector<Group> groups;
            std::vector<bool> assigned(ops.size(), false);

            for (size_t i = 0; i < ops.size(); ++i)
            {
                if (assigned[i] || ops[i].get_type() != Operation::Type::CCX)
                    continue;

                const auto &q = ops[i].get_qubits();
                if (q.size() != 3 || q[0] == q[1] || q[0] == q[2] || q[1] == q[2])
                {
                    throw std::runtime_error(
                        "CCXLoweringPass: CCX requires three distinct qubit operands");
                }
                if (ops[i].is_conditional())
                {
                    throw std::runtime_error(
                        "CCXLoweringPass: lowering a classically conditioned CCX is not "
                        "supported; the gadget's own measurement would need conditioning too");
                }

                Group g;
                g.c0 = q[0];
                g.c1 = q[1];
                g.members.push_back(i);
                assigned[i] = true;

                for (size_t j = i + 1; j < ops.size(); ++j)
                {
                    const auto &oj = ops[j];
                    if (!assigned[j] && oj.get_type() == Operation::Type::CCX &&
                        !oj.is_conditional() && oj.get_qubits().size() == 3 &&
                        same_controls(oj, g))
                    {
                        g.members.push_back(j);
                        assigned[j] = true;
                        continue;
                    }
                    auto hit = disturbed_qubits(oj);
                    if (std::find(hit.begin(), hit.end(), g.c0) != hit.end() ||
                        std::find(hit.begin(), hit.end(), g.c1) != hit.end())
                    {
                        break;
                    }
                }
                groups.push_back(std::move(g));
            }
            return groups;
        }

        static bool same_controls(const Operation &op, const Group &g)
        {
            const auto &q = op.get_qubits();
            return (q[0] == g.c0 && q[1] == g.c1) || (q[0] == g.c1 && q[1] == g.c0);
        }

        static void split_to_singletons(std::vector<Group> &groups)
        {
            std::vector<Group> out;
            out.reserve(groups.size());
            for (const auto &g : groups)
            {
                for (size_t m : g.members)
                {
                    Group s;
                    s.c0 = g.c0;
                    s.c1 = g.c1;
                    s.members = {m};
                    out.push_back(std::move(s));
                }
            }
            std::sort(out.begin(), out.end(),
                      [](const Group &a, const Group &b) { return a.begin() < b.begin(); });
            groups = std::move(out);
        }

        /** Halve the widest group until no more than `budget` are ever live at once. */
        static void enforce_budget(std::vector<Group> &groups, size_t budget)
        {
            while (max_overlap(groups) > budget)
            {
                auto widest = std::max_element(
                    groups.begin(), groups.end(), [](const Group &a, const Group &b) {
                        return a.members.size() < b.members.size();
                    });
                if (widest == groups.end() || widest->members.size() < 2)
                    break; // all singletons already; cannot reduce further
                Group tail;
                tail.c0 = widest->c0;
                tail.c1 = widest->c1;
                const size_t half = widest->members.size() / 2;
                tail.members.assign(widest->members.begin() + half, widest->members.end());
                widest->members.resize(half);
                groups.push_back(std::move(tail));
                std::sort(groups.begin(), groups.end(),
                          [](const Group &a, const Group &b) { return a.begin() < b.begin(); });
            }
        }

        static size_t max_overlap(const std::vector<Group> &groups)
        {
            std::vector<std::pair<size_t, int>> ev;
            ev.reserve(groups.size() * 2);
            for (const auto &g : groups)
            {
                ev.emplace_back(g.begin(), +1);
                ev.emplace_back(g.end(), -1);
            }
            // Open before close at equal positions, so a single-member group whose
            // begin equals its end still counts as live for one step. Distinct groups
            // never share a position, since each operation index belongs to at most one
            // member, so this cannot create a false overlap.
            std::sort(ev.begin(), ev.end(), [](const auto &a, const auto &b) {
                return a.first != b.first ? a.first < b.first : a.second > b.second;
            });
            size_t live = 0, mx = 0;
            for (const auto &e : ev)
            {
                if (e.second > 0)
                    mx = std::max(mx, ++live);
                else
                    --live;
            }
            return mx;
        }

        /**
         * @brief Interval-graph colouring: a greedy sweep by start position uses exactly
         *        the minimum number of ancillas, which is the maximum group overlap.
         */
        static size_t assign_ancillas(std::vector<Group> &groups)
        {
            std::sort(groups.begin(), groups.end(),
                      [](const Group &a, const Group &b) { return a.begin() < b.begin(); });
            std::vector<size_t> free_at; // free_at[c] = first position colour c is free
            for (auto &g : groups)
            {
                size_t chosen = free_at.size();
                for (size_t c = 0; c < free_at.size(); ++c)
                {
                    if (free_at[c] <= g.begin())
                    {
                        chosen = c;
                        break;
                    }
                }
                if (chosen == free_at.size())
                    free_at.push_back(0);
                free_at[chosen] = g.end() + 1;
                g.ancilla = chosen;
            }
            return free_at.size();
        }

        void emit_and_compute(const Group &g, Circuit &out) const
        {
            for (auto &op : and_compute(g.c0, g.c1, g.ancilla))
                out.add_operation(std::move(op));
        }

        void emit_and_uncompute(const Group &g, Circuit &out) const
        {
            for (auto &op : and_uncompute(g.c0, g.c1, g.ancilla, g.result_bit, options_.emit_reset))
                out.add_operation(std::move(op));
        }

        CCXLoweringOptions options_;
        std::vector<Group> groups_;
        size_t ancillas_added_ = 0;
        size_t result_bits_added_ = 0;
    };

} // namespace NWQEC
