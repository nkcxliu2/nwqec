#pragma once

#include "pass_template.hpp"
#include "and_gadget.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace NWQEC
{
    /**
     * @brief How a multi-controlled X is realized.
     *
     * Every strategy here uses only *clean* ancillae: each enters in |0> and is returned
     * to |0>. They span a space-time trade-off rather than a ranking, and none is the
     * default -- selection is left to the caller, and eventually to a combined
     * space-time cost model.
     *
     * With n controls, and reading T-counts under the downstream CCXStrategy::Standard7T:
     *
     *   AndLadderClean   n-1 ancillas, 4n-4 T,   O(n) depth
     *   AndTreeClean     n-1 ancillas, 4n-4 T,   O(log n) depth
     *   VChainM15        n-2 ancillas, 4n-1 T,   O(n) depth
     *   HybridBudget     A   ancillas, 8n-4A-5 T (A free), O(n) depth
     *   KG24_1Clean      1   ancilla,  8n-9 T,   O(n) depth
     *   KG24_2Clean      2   ancillas, 8n-9 T,   O(log n) depth
     */
    enum class MCXStrategy
    {
        Preserve,       // leave MCX untouched
        AndLadderClean, // linear chain of Gidney ANDs
        AndTreeClean,   // balanced tree of Gidney ANDs
        VChainM15,      // Maslov V-chain: AND chain with a real Toffoli at the target
        HybridBudget,   // AND prefix within an ancilla budget, KG24 for the remainder
        KG24_1Clean,    // Khattar-Gidney conditionally clean ancillae, linear depth
        KG24_2Clean     // Khattar-Gidney conditionally clean ancillae, log depth
    };

    struct MCXLoweringOptions
    {
        MCXStrategy strategy = MCXStrategy::Preserve;
        // Ancilla budget, and the budget A for HybridBudget. Following the convention of
        // CCXLoweringPass, 0 means unlimited -- which for HybridBudget means A = n-1 and
        // hence degeneration to AndLadderClean.
        size_t max_ancillas = 0;
        // Whether an ancilla may be used at all. No clean-ancilla strategy applies when
        // this is false, so an MCX with three or more controls is left in place.
        bool allow_ancilla = true;
        bool emit_reset = true;
    };

    /** Resource footprint of one lowering, in gate classes rather than T directly. */
    struct MCXCost
    {
        size_t and_count = 0;  // Gidney ANDs: 4 T each, erased for free
        size_t rccx_count = 0; // bare RCCX in KG24 ladders: 4 T each, no erasure available
        size_t ccx_count = 0;  // terminal Toffolis, handed to CCXLoweringPass
        size_t clean_ancillas = 0;
        size_t result_bits = 0;
        const char *depth_class = "O(n)";

        /** T under CCXStrategy::Standard7T (7 T per terminal Toffoli, no extra ancilla). */
        size_t t_standard7t() const { return 4 * and_count + 4 * rccx_count + 7 * ccx_count; }
        /** T under CCXStrategy::Jones4T (4 T per terminal Toffoli, one ancilla each). */
        size_t t_jones4t() const { return 4 * and_count + 4 * rccx_count + 4 * ccx_count; }
        size_t ancillas_jones4t() const { return clean_ancillas + ccx_count; }
    };

    /**
     * @brief Lowers MCX to RCCX / CCX / CX / X using clean ancillae.
     *
     * Must run before CCX_LOWERING: the terminal Toffolis emitted by VChainM15,
     * HybridBudget and the KG24 strategies are handed on for that pass to realize, so
     * the T-count of every strategy depends on the downstream CCXStrategy.
     *
     * Ancillae are live only within the span of a single MCX, so one pool sized to the
     * largest single requirement is reused across every MCX in the circuit.
     *
     * References:
     *   C. Gidney, Quantum 2, 74 (2018)                     -- temporary logical AND
     *   D. Maslov, Phys. Rev. A 93, 022311 (2016)           -- V-chain
     *   N. Khattar and C. Gidney, Quantum 9, 1752 (2025)    -- conditionally clean ancillae
     */
    class MCXLoweringPass : public Pass
    {
    public:
        explicit MCXLoweringPass(MCXLoweringOptions options = {}) : options_(options) {}

        std::string get_name() const override { return "MCX Lowering Pass"; }

        size_t ancillas_added() const { return ancillas_added_; }
        size_t result_bits_added() const { return result_bits_added_; }

        // ---- static cost model -------------------------------------------------
        //
        // Derived by building the expansion and counting it, so the model cannot drift
        // away from what the pass actually emits.

        static MCXCost mcx_cost(size_t num_controls, MCXStrategy strategy, size_t budget = 0)
        {
            MCXCost cost;
            if (strategy == MCXStrategy::Preserve)
                return cost;
            cost.depth_class = (strategy == MCXStrategy::AndTreeClean ||
                                strategy == MCXStrategy::KG24_2Clean)
                                   ? "O(log n)"
                                   : "O(n)";

            const size_t n_anc = ancilla_requirement(num_controls, strategy, budget);
            std::vector<size_t> anc(n_anc), bits(n_anc);
            // Indices past the controls and target, so they cannot collide with them.
            std::iota(anc.begin(), anc.end(), num_controls + 1);
            std::iota(bits.begin(), bits.end(), 0);

            std::vector<size_t> controls(num_controls);
            std::iota(controls.begin(), controls.end(), 0);

            std::vector<Operation> ops;
            expand(controls, num_controls, anc, bits, strategy, budget, true, ops);

            std::vector<bool> anc_used(n_anc, false);
            for (const auto &op : ops)
            {
                switch (op.get_type())
                {
                case Operation::Type::RCCX:
                    // An AND is RCCX followed by SDG; the SDG is what distinguishes it.
                    ++cost.rccx_count;
                    break;
                case Operation::Type::SDG:
                    --cost.rccx_count;
                    ++cost.and_count;
                    break;
                case Operation::Type::CCX:
                    ++cost.ccx_count;
                    break;
                case Operation::Type::MEASURE:
                    ++cost.result_bits;
                    break;
                default:
                    break;
                }
                for (size_t q : op.get_qubits())
                {
                    if (q > num_controls && q - num_controls - 1 < n_anc)
                        anc_used[q - num_controls - 1] = true;
                }
            }
            cost.clean_ancillas =
                static_cast<size_t>(std::count(anc_used.begin(), anc_used.end(), true));
            cost.result_bits = std::min(cost.result_bits, cost.clean_ancillas);
            return cost;
        }

        /** Ancillae the strategy needs for n controls; the pool is sized by the max. */
        static size_t ancilla_requirement(size_t n, MCXStrategy strategy, size_t budget)
        {
            if (n < 3 || strategy == MCXStrategy::Preserve)
                return 0;
            switch (strategy)
            {
            case MCXStrategy::AndLadderClean:
            case MCXStrategy::AndTreeClean:
                return n - 1;
            case MCXStrategy::VChainM15:
                return n - 2;
            case MCXStrategy::KG24_1Clean:
                return 1;
            case MCXStrategy::KG24_2Clean:
                return kg24_2_ancillas(n);
            case MCXStrategy::HybridBudget:
                return hybrid_budget(n, budget);
            default:
                return 0;
            }
        }

        /** Bits the strategy needs; only the AND-based ones measure. */
        static size_t bit_requirement(size_t n, MCXStrategy strategy, size_t budget)
        {
            if (n < 3 || strategy == MCXStrategy::Preserve)
                return 0;
            switch (strategy)
            {
            case MCXStrategy::AndLadderClean:
            case MCXStrategy::AndTreeClean:
                return n - 1;
            case MCXStrategy::VChainM15:
                return n - 2;
            case MCXStrategy::HybridBudget:
            {
                const size_t a = hybrid_budget(n, budget);
                return (a >= n - 1) ? n - 1 : a - 1;
            }
            default:
                return 0;
            }
        }

        bool run(Circuit &circuit) override
        {
            if (options_.strategy == MCXStrategy::Preserve)
                return false;

            const auto &ops = circuit.get_operations();
            size_t need_anc = 0, need_bits = 0;
            bool any_mcx = false, any_needs_ancilla = false;
            for (const auto &op : ops)
            {
                if (op.get_type() != Operation::Type::MCX)
                    continue;
                validate(op);
                any_mcx = true;
                const size_t n = op.get_qubits().size() - 1;
                if (n < 3)
                    continue; // CX and CCX need no ancilla
                any_needs_ancilla = true;
                const MCXStrategy s = effective_strategy(n);
                if (s == MCXStrategy::Preserve)
                    continue;
                need_anc = std::max(need_anc, ancilla_requirement(n, s, options_.max_ancillas));
                need_bits = std::max(need_bits, bit_requirement(n, s, options_.max_ancillas));
            }
            if (!any_mcx)
                return false;

            if (!options_.allow_ancilla && any_needs_ancilla)
            {
                std::cerr << "MCXLoweringPass: allow_ancilla is false and no ancilla-free MCX "
                             "construction is available; leaving MCX with three or more controls "
                             "in place. Downstream passes will not understand them.\n";
                need_anc = 0;
                need_bits = 0;
            }

            const size_t nq_in = circuit.get_num_qubits();
            const size_t nb_in = circuit.get_num_bits();

            Circuit out;
            out.add_qreg("q", nq_in + need_anc);
            if (nb_in + need_bits > 0)
                out.add_creg("c", nb_in + need_bits);

            std::vector<size_t> anc(need_anc), bits(need_bits);
            std::iota(anc.begin(), anc.end(), nq_in);
            std::iota(bits.begin(), bits.end(), nb_in);

            bool modified = false;
            for (const auto &op : ops)
            {
                if (op.get_type() != Operation::Type::MCX)
                {
                    out.add_operation(op);
                    continue;
                }
                const auto &q = op.get_qubits();
                const size_t n = q.size() - 1;
                const std::vector<size_t> controls(q.begin(), q.end() - 1);
                const size_t target = q.back();

                MCXStrategy s = (n < 3) ? options_.strategy : effective_strategy(n);
                if (n >= 3 && (!options_.allow_ancilla || s == MCXStrategy::Preserve))
                {
                    out.add_operation(op); // could not be lowered; diagnostic already issued
                    continue;
                }

                std::vector<Operation> body;
                expand(controls, target, anc, bits, s, options_.max_ancillas,
                       options_.emit_reset, body);
                apply_condition(op, body);
                for (auto &b : body)
                    out.add_operation(std::move(b));
                modified = true;
            }

            if (!modified)
                return false;

            ancillas_added_ = need_anc;
            result_bits_added_ = need_bits;
            circuit = std::move(out);
            return true;
        }

        // ---- expansion ---------------------------------------------------------

        /**
         * @brief Emit the chosen realization of C^nX into @p out.
         *
         * @p anc and @p bits supply the ancilla qubits and their measurement bits; they
         * must be at least as long as ancilla_requirement / bit_requirement report.
         */
        static void expand(const std::vector<size_t> &controls, size_t target,
                           const std::vector<size_t> &anc, const std::vector<size_t> &bits,
                           MCXStrategy strategy, size_t budget, bool emit_reset,
                           std::vector<Operation> &out)
        {
            const size_t n = controls.size();
            if (n == 0)
            {
                out.push_back(Operation(Operation::Type::X, {target}));
                return;
            }
            if (n == 1)
            {
                out.push_back(Operation(Operation::Type::CX, {controls[0], target}));
                return;
            }
            if (n == 2)
            {
                out.push_back(Operation(Operation::Type::CCX, {controls[0], controls[1], target}));
                return;
            }

            switch (strategy)
            {
            case MCXStrategy::AndLadderClean:
                emit_and_ladder(controls, target, anc, bits, emit_reset, out);
                break;
            case MCXStrategy::AndTreeClean:
                emit_and_tree(controls, target, anc, bits, emit_reset, out);
                break;
            case MCXStrategy::VChainM15:
                emit_vchain(controls, target, anc, bits, emit_reset, out);
                break;
            case MCXStrategy::HybridBudget:
                emit_hybrid(controls, target, anc, bits, budget, emit_reset, out);
                break;
            case MCXStrategy::KG24_1Clean:
                emit_kg24_1(controls, target, anc[0], out);
                break;
            case MCXStrategy::KG24_2Clean:
                emit_kg24_2(controls, target, anc, out);
                break;
            default:
                throw std::runtime_error("MCXLoweringPass: no expansion for this strategy");
            }
        }

    private:
        struct AndRec
        {
            size_t a, b, w, bit;
        };

        static void append(std::vector<Operation> &dst, std::vector<Operation> src)
        {
            for (auto &op : src)
                dst.push_back(std::move(op));
        }

        static void unwind(const std::vector<AndRec> &recs, bool emit_reset,
                           std::vector<Operation> &out)
        {
            for (auto it = recs.rbegin(); it != recs.rend(); ++it)
                append(out, and_uncompute(it->a, it->b, it->w, it->bit, emit_reset));
        }

        /**
         * @brief Fold the first @p take controls into one ancilla with a chain of ANDs.
         *
         * w_0 = c_0 AND c_1, w_i = w_{i-1} AND c_{i+1}, consuming take-1 ancillae and
         * leaving the conjunction of controls[0..take-1] in the returned qubit.
         */
        static size_t emit_prefix_ladder(const std::vector<size_t> &controls, size_t take,
                                         const std::vector<size_t> &anc,
                                         const std::vector<size_t> &bits,
                                         std::vector<AndRec> &recs, std::vector<Operation> &out)
        {
            size_t head = anc[0];
            recs.push_back({controls[0], controls[1], anc[0], bits[0]});
            append(out, and_compute(controls[0], controls[1], anc[0]));
            for (size_t i = 1; i + 1 < take; ++i)
            {
                recs.push_back({head, controls[i + 1], anc[i], bits[i]});
                append(out, and_compute(head, controls[i + 1], anc[i]));
                head = anc[i];
            }
            return head;
        }

        static void emit_and_ladder(const std::vector<size_t> &controls, size_t target,
                                    const std::vector<size_t> &anc, const std::vector<size_t> &bits,
                                    bool emit_reset, std::vector<Operation> &out)
        {
            std::vector<AndRec> recs;
            const size_t head = emit_prefix_ladder(controls, controls.size(), anc, bits, recs, out);
            out.push_back(Operation(Operation::Type::CX, {head, target}));
            unwind(recs, emit_reset, out);
        }

        static void emit_vchain(const std::vector<size_t> &controls, size_t target,
                                const std::vector<size_t> &anc, const std::vector<size_t> &bits,
                                bool emit_reset, std::vector<Operation> &out)
        {
            std::vector<AndRec> recs;
            const size_t n = controls.size();
            const size_t head = emit_prefix_ladder(controls, n - 1, anc, bits, recs, out);
            out.push_back(Operation(Operation::Type::CCX, {head, controls[n - 1], target}));
            unwind(recs, emit_reset, out);
        }

        static void emit_and_tree(const std::vector<size_t> &controls, size_t target,
                                  const std::vector<size_t> &anc, const std::vector<size_t> &bits,
                                  bool emit_reset, std::vector<Operation> &out)
        {
            std::vector<AndRec> recs;
            std::vector<size_t> level = controls;
            size_t next = 0;
            while (level.size() > 1)
            {
                std::vector<size_t> up;
                size_t i = 0;
                for (; i + 1 < level.size(); i += 2)
                {
                    const size_t w = anc[next];
                    recs.push_back({level[i], level[i + 1], w, bits[next]});
                    append(out, and_compute(level[i], level[i + 1], w));
                    up.push_back(w);
                    ++next;
                }
                if (i < level.size())
                    up.push_back(level[i]); // odd one out, carried to the next layer
                level = std::move(up);
            }
            out.push_back(Operation(Operation::Type::CX, {level[0], target}));
            unwind(recs, emit_reset, out);
        }

        static void emit_hybrid(const std::vector<size_t> &controls, size_t target,
                                const std::vector<size_t> &anc, const std::vector<size_t> &bits,
                                size_t budget, bool emit_reset, std::vector<Operation> &out)
        {
            const size_t n = controls.size();
            const size_t a = hybrid_budget(n, budget);
            if (a >= n - 1)
            {
                // At full budget the pure AND ladder is strictly better than folding a
                // prefix and finishing with KG24: it ends in an AND (4 T) rather than a
                // terminal Toffoli (7 T).
                emit_and_ladder(controls, target, anc, bits, emit_reset, out);
                return;
            }

            std::vector<AndRec> recs;
            std::vector<size_t> rest;
            if (a == 1)
            {
                rest = controls; // nothing folded; the whole problem goes to KG24
            }
            else
            {
                const size_t head = emit_prefix_ladder(controls, a, anc, bits, recs, out);
                rest.push_back(head);
                rest.insert(rest.end(), controls.begin() + static_cast<long>(a), controls.end());
            }
            emit_kg24_1(rest, target, anc[a - 1], out);
            unwind(recs, emit_reset, out);
        }

        /**
         * @brief Khattar-Gidney ladder, Fig. 3 of arXiv:2407.17966.
         *
         * @p map takes local indices to circuit qubits: local 0 is the ancilla and local
         * 1+j is control j. Transcribed against Qiskit's synth_mcx_1_clean_kg24.
         */
        static void emit_kg24_ladder(long long N, const std::vector<size_t> &map,
                                     std::vector<Operation> &out)
        {
            auto rccx = [&](long long x, long long y, long long t) {
                out.push_back(Operation(Operation::Type::RCCX,
                                        {map[static_cast<size_t>(x)], map[static_cast<size_t>(y)],
                                         map[static_cast<size_t>(t)]}));
                out.push_back(Operation(Operation::Type::X, {map[static_cast<size_t>(t)]}));
            };

            for (long long i = 2; i < N - 2; i += 2)
                rccx(i + 1, i + 2, i);

            long long a, b, tgt;
            if (N % 2 != 0)
            {
                a = N - 3;
                b = N - 5;
                tgt = N - 6;
            }
            else
            {
                a = N - 1;
                b = N - 4;
                tgt = N - 5;
            }
            if (tgt > 0)
                rccx(a, b, tgt);
            for (long long i = tgt; i > 2; i -= 2)
                rccx(i, i - 1, i - 2);
        }

        static void emit_kg24_1(const std::vector<size_t> &controls, size_t target, size_t anc,
                                std::vector<Operation> &out)
        {
            const size_t n = controls.size();
            if (n < 3)
            {
                out.push_back(Operation(Operation::Type::CCX, {controls[0], controls[1], target}));
                return;
            }
            std::vector<size_t> map;
            map.reserve(n + 1);
            map.push_back(anc);
            map.insert(map.end(), controls.begin(), controls.end());

            const long long N = static_cast<long long>(n) + 1;
            std::vector<Operation> ladder;
            emit_kg24_ladder(N, map, ladder);
            const size_t final_ctrl = static_cast<size_t>(std::max(0LL, 6LL - N));

            out.push_back(Operation(Operation::Type::RCCX, {controls[0], controls[1], anc}));
            out.insert(out.end(), ladder.begin(), ladder.end());
            out.push_back(Operation(Operation::Type::CCX, {anc, controls[final_ctrl], target}));
            // RCCX and X are self-inverse in NWQEC, so the inverse ladder is the reverse.
            out.insert(out.end(), ladder.rbegin(), ladder.rend());
            out.push_back(Operation(Operation::Type::RCCX, {controls[0], controls[1], anc}));
        }

        /**
         * @brief Log-depth conditionally clean ladder, Fig. 4b of arXiv:2407.17966.
         *
         * Works in local indices 0..n-1 for controls and n for the ancilla; returns the
         * controls left over for the linear-depth call. Transcribed against Qiskit's
         * _build_logn_depth_ccx_ladder.
         */
        static void emit_logn_ladder(size_t n, const std::vector<size_t> &map,
                                     std::vector<Operation> &out,
                                     std::vector<size_t> &final_ctrls_out)
        {
            const size_t anc_local = n;
            std::vector<size_t> ctrls(n);
            std::iota(ctrls.begin(), ctrls.end(), 0);
            std::vector<size_t> anc{anc_local};
            std::vector<size_t> final_ctrls;

            while (ctrls.size() > 1)
            {
                const size_t take = std::min(anc.size() + 1, ctrls.size());
                std::vector<size_t> batch(ctrls.begin(), ctrls.begin() + static_cast<long>(take));
                ctrls.erase(ctrls.begin(), ctrls.begin() + static_cast<long>(take));

                std::vector<size_t> new_anc;
                while (batch.size() > 1)
                {
                    const size_t m = batch.size() / 2;
                    const size_t st = batch.size() % 2;
                    const std::vector<size_t> xs(batch.begin() + static_cast<long>(st),
                                                 batch.begin() + static_cast<long>(st + m));
                    const std::vector<size_t> ys(batch.begin() + static_cast<long>(st + m),
                                                 batch.end());
                    const std::vector<size_t> ts(anc.end() - static_cast<long>(m), anc.end());

                    if (ts.size() == 1 && ts[0] == anc_local)
                    {
                        // The one genuinely clean ancilla; no X, matching Fig. 4b step 1.
                        out.push_back(Operation(Operation::Type::RCCX,
                                                {map[xs[0]], map[ys[0]], map[ts[0]]}));
                    }
                    else
                    {
                        // Here the X *precedes* the RCCX, unlike the linear ladder.
                        for (size_t k = 0; k < m; ++k)
                            out.push_back(Operation(Operation::Type::X, {map[ts[k]]}));
                        for (size_t k = 0; k < m; ++k)
                            out.push_back(Operation(Operation::Type::RCCX,
                                                    {map[xs[k]], map[ys[k]], map[ts[k]]}));
                    }

                    new_anc.insert(new_anc.end(), batch.begin() + static_cast<long>(st),
                                   batch.end());
                    std::vector<size_t> nb = ts;
                    nb.insert(nb.end(), batch.begin(), batch.begin() + static_cast<long>(st));
                    batch = std::move(nb);
                    anc.resize(anc.size() - m);
                }
                anc.insert(anc.end(), new_anc.begin(), new_anc.end());
                std::sort(anc.begin(), anc.end());
                final_ctrls.insert(final_ctrls.end(), batch.begin(), batch.end());
            }
            final_ctrls.insert(final_ctrls.end(), ctrls.begin(), ctrls.end());
            std::sort(final_ctrls.begin(), final_ctrls.end());
            final_ctrls.pop_back(); // the ancilla, which sorts last
            final_ctrls_out = std::move(final_ctrls);
        }

        static void emit_kg24_2(const std::vector<size_t> &controls, size_t target,
                                const std::vector<size_t> &anc, std::vector<Operation> &out)
        {
            const size_t n = controls.size();
            std::vector<size_t> map(controls);
            map.push_back(anc[0]);

            std::vector<Operation> ladder;
            std::vector<size_t> fc;
            emit_logn_ladder(n, map, ladder, fc);

            out.insert(out.end(), ladder.begin(), ladder.end());
            if (fc.size() == 1)
            {
                out.push_back(Operation(Operation::Type::CCX, {anc[0], controls[fc[0]], target}));
            }
            else
            {
                std::vector<size_t> mid{anc[0]};
                for (size_t i : fc)
                    mid.push_back(controls[i]);
                emit_kg24_1(mid, target, anc[1], out);
            }
            out.insert(out.end(), ladder.rbegin(), ladder.rend());
        }

        // ---- helpers -----------------------------------------------------------

        static size_t hybrid_budget(size_t n, size_t budget)
        {
            if (budget == 0)
                return n - 1; // unlimited, per the CCXLoweringPass convention
            return std::max<size_t>(1, std::min(budget, n - 1));
        }

        static size_t kg24_2_ancillas(size_t n)
        {
            // The second ancilla is only needed when the log-depth ladder leaves more
            // than one control behind; for small n it leaves exactly one.
            std::vector<size_t> map(n + 1);
            std::iota(map.begin(), map.end(), 0);
            std::vector<Operation> scratch;
            std::vector<size_t> fc;
            emit_logn_ladder(n, map, scratch, fc);
            return fc.size() == 1 ? 1 : 2;
        }

        /** Fall back down the trade-off table to the cheapest-in-space strategy that fits. */
        MCXStrategy effective_strategy(size_t n) const
        {
            const size_t budget = options_.max_ancillas;
            if (budget == 0 || options_.strategy == MCXStrategy::HybridBudget)
                return options_.strategy; // unlimited, or a strategy that clamps itself
            if (ancilla_requirement(n, options_.strategy, budget) <= budget)
                return options_.strategy;
            for (MCXStrategy f : {MCXStrategy::VChainM15, MCXStrategy::KG24_2Clean,
                                  MCXStrategy::KG24_1Clean})
            {
                if (ancilla_requirement(n, f, budget) <= budget)
                    return f;
            }
            return MCXStrategy::Preserve;
        }

        static void validate(const Operation &op)
        {
            const auto &q = op.get_qubits();
            if (q.size() < 2)
                throw std::runtime_error("MCXLoweringPass: MCX needs at least one control "
                                         "and a target");
            std::vector<size_t> sorted(q);
            std::sort(sorted.begin(), sorted.end());
            if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
                throw std::runtime_error("MCXLoweringPass: MCX operands must be distinct");
        }

        /**
         * @brief Propagate a classical guard onto the expansion.
         *
         * The KG24 strategies emit no measurement, so guarding them is just guarding
         * every gate. The AND-based ones measure, and a measurement cannot be guarded,
         * so a conditional MCX is rejected there rather than silently mis-lowered.
         */
        static void apply_condition(const Operation &src, std::vector<Operation> &body)
        {
            if (!src.is_conditional())
                return;
            for (const auto &op : body)
            {
                if (!Operation::is_conditionable(op.get_type()))
                {
                    throw std::runtime_error(
                        "MCXLoweringPass: cannot lower a classically conditioned MCX with a "
                        "measurement-based strategy; use KG24_1Clean or KG24_2Clean, which "
                        "emit no measurement");
                }
            }
            for (auto &op : body)
                op.set_condition(*src.get_condition());
        }

        MCXLoweringOptions options_;
        size_t ancillas_added_ = 0;
        size_t result_bits_added_ = 0;
    };

} // namespace NWQEC
