// Python bindings for NWQEC using pybind11

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

#include "nwqec/parser/qasm_parser.hpp"
#include "nwqec/core/operation.hpp"
#include "nwqec/core/pauli_op.hpp"
#include "nwqec/core/constants.hpp"

#include "nwqec/core/transpiler.hpp"
#include "nwqec/analysis/clifford_t_counts.hpp"

namespace py = pybind11;

namespace
{
    // Lightweight, non-owning handle used by Circuit.c_if(...). It appends gates through
    // the circuit's ordinary path, then guards whatever was appended, so a conditional
    // gate follows exactly the same construction rules as an unconditional one.
    struct ConditionalBuilder
    {
        NWQEC::Circuit *circuit;
        NWQEC::ClassicalCondition condition;

        template <typename F>
        ConditionalBuilder &guarded(F &&emit)
        {
            size_t first = circuit->get_operations().size();
            emit(*circuit);
            circuit->set_condition_on_operations(first, condition);
            return *this;
        }
    };

    // Helper to render stats to a string
    std::string circuit_stats(const NWQEC::Circuit &c)
    {
        std::ostringstream oss;
        c.print_stats(oss);
        return oss.str();
    }
    std::string circuit_to_qasm(const NWQEC::Circuit &c)
    {
        std::ostringstream oss;
        c.print(oss);
        return oss.str();
    }
    void circuit_save_qasm(const NWQEC::Circuit &c, const std::string &filename)
    {
        std::ofstream ofs(filename);
        if (!ofs)
        {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        c.print(ofs);
    }

    // Helpers to enforce PBC vs standard gate exclusivity in Python API
    inline bool is_pauli_op(NWQEC::Operation::Type t)
    {
        using T = NWQEC::Operation::Type;
        return t == T::T_PAULI || t == T::S_PAULI || t == T::Z_PAULI || t == T::M_PAULI;
    }

    inline bool is_barrier(NWQEC::Operation::Type t)
    {
        return t == NWQEC::Operation::Type::BARRIER;
    }

    bool circuit_has_pauli_ops(const NWQEC::Circuit &c)
    {
        for (const auto &op : c.get_operations())
        {
            if (is_pauli_op(op.get_type()))
                return true;
        }
        return false;
    }

    bool circuit_has_non_pauli_ops(const NWQEC::Circuit &c)
    {
        for (const auto &op : c.get_operations())
        {
            if (!is_pauli_op(op.get_type()) && !is_barrier(op.get_type()))
                return true;
        }
        return false;
    }

    NWQEC::RzErrorPolicy parse_rz_err_policy(const std::string &rz_err)
    {
        if (rz_err == "per-gate")
            return NWQEC::RzErrorPolicy::PER_GATE;
        if (rz_err == "total")
            return NWQEC::RzErrorPolicy::TOTAL;
        if (rz_err == "relative")
            return NWQEC::RzErrorPolicy::RELATIVE;

        throw std::invalid_argument("rz_err must be one of: per-gate, total, relative");
    }

    NWQEC::CCXStrategy parse_ccx_strategy(const std::string &name)
    {
        if (name == "standard-7t")
            return NWQEC::CCXStrategy::Standard7T;
        if (name == "jones-4t")
            return NWQEC::CCXStrategy::Jones4T;
        if (name == "shared-and")
            return NWQEC::CCXStrategy::SharedAND;
        if (name == "preserve")
            return NWQEC::CCXStrategy::Preserve;
        throw std::invalid_argument(
            "ccx_strategy must be one of: standard-7t, jones-4t, shared-and, preserve");
    }

    // Resolve the CCX strategy, rejecting an inconsistent combination with the
    // legacy keep_ccx flag rather than silently picking one.
    // An empty strategy string means "not specified", so the legacy keep_ccx flag still
    // decides on its own. Only an explicitly given strategy can conflict with it.
    NWQEC::CCXLoweringOptions resolve_ccx_options(const std::string &strategy,
                                                  bool keep_ccx,
                                                  size_t max_ancillas,
                                                  bool allow_ancilla)
    {
        NWQEC::CCXLoweringOptions opts;
        opts.max_ancillas = max_ancillas;
        opts.allow_ancilla = allow_ancilla;

        if (strategy.empty())
        {
            opts.strategy = keep_ccx ? NWQEC::CCXStrategy::Preserve
                                     : NWQEC::CCXStrategy::Standard7T;
            return opts;
        }

        opts.strategy = parse_ccx_strategy(strategy);
        if (keep_ccx && opts.strategy != NWQEC::CCXStrategy::Preserve)
        {
            throw std::invalid_argument(
                "keep_ccx=True conflicts with ccx_strategy='" + strategy +
                "'; use ccx_strategy='preserve' instead");
        }
        return opts;
    }

    NWQEC::MCXStrategy parse_mcx_strategy(const std::string &name)
    {
        if (name == "preserve")
            return NWQEC::MCXStrategy::Preserve;
        if (name == "and-ladder")
            return NWQEC::MCXStrategy::AndLadderClean;
        if (name == "and-tree")
            return NWQEC::MCXStrategy::AndTreeClean;
        if (name == "vchain-m15")
            return NWQEC::MCXStrategy::VChainM15;
        if (name == "hybrid-budget")
            return NWQEC::MCXStrategy::HybridBudget;
        if (name == "kg24-1-clean")
            return NWQEC::MCXStrategy::KG24_1Clean;
        if (name == "kg24-2-clean")
            return NWQEC::MCXStrategy::KG24_2Clean;
        throw std::invalid_argument(
            "mcx_strategy must be one of: preserve, and-ladder, and-tree, vchain-m15, "
            "hybrid-budget, kg24-1-clean, kg24-2-clean");
    }

    NWQEC::MCXLoweringOptions resolve_mcx_options(const std::string &strategy,
                                                  size_t max_ancillas,
                                                  bool allow_ancilla)
    {
        NWQEC::MCXLoweringOptions opts;
        opts.max_ancillas = max_ancillas;
        opts.allow_ancilla = allow_ancilla;
        opts.strategy = strategy.empty() ? NWQEC::MCXStrategy::Preserve
                                         : parse_mcx_strategy(strategy);
        return opts;
    }

    void configure_rz_error(NWQEC::PassConfig &config, const std::string &rz_err, py::object epsilon)
    {
        config.rz_error_policy = parse_rz_err_policy(rz_err);
        double resolved_epsilon = epsilon.is_none()
                                      ? NWQEC::default_rz_error_epsilon(config.rz_error_policy)
                                      : epsilon.cast<double>();

        if (!(resolved_epsilon > 0.0) || !std::isfinite(resolved_epsilon))
        {
            throw std::invalid_argument("epsilon must be a positive finite number");
        }

        config.rz_error_epsilon = resolved_epsilon;
    }

    // Helper to run transforms using the Transpiler
    std::unique_ptr<NWQEC::Circuit> apply_transforms(const NWQEC::Circuit &circuit,
                                                     bool to_pbc,
                                                     bool to_clifford_reduction,
                                                     bool keep_cx,
                                                     bool t_pauli_opt,
                                                     bool remove_pauli,
                                                     bool keep_ccx,
                                                     bool silent,
                                                     const std::string &rz_err,
                                                     py::object epsilon,
                                                     NWQEC::CCXLoweringOptions ccx_options = {},
                                                     NWQEC::MCXLoweringOptions mcx_options = {})
    {
        NWQEC::Transpiler transpiler;
        NWQEC::PassConfig config;
        // DecomposePass reads the boolean, so the Preserve strategy has to reach it too.
        config.keep_ccx = keep_ccx || (ccx_options.strategy == NWQEC::CCXStrategy::Preserve);
        config.ccx_lowering = ccx_options;
        config.mcx_lowering = mcx_options;
        config.keep_cx = keep_cx;
        config.silent = silent;
        configure_rz_error(config, rz_err, epsilon);
        
        auto circuit_copy = std::make_unique<NWQEC::Circuit>(circuit);
        
        // Choose the appropriate pass sequence
        std::vector<NWQEC::PassType> passes;
        
        if (t_pauli_opt) {
            // T-optimization only - assumes input is already PBC
            passes = NWQEC::PassSequences::T_OPTIMIZATION_ONLY;
        } else if (to_pbc) {
            passes = NWQEC::PassSequences::TO_PBC;
        } else if (to_clifford_reduction) {
            passes = NWQEC::PassSequences::TO_CLIFFORD_REDUCTION;
        } else {
            // Default: Clifford+T conversion
            passes = NWQEC::PassSequences::TO_CLIFFORD_T;
        }
        
        // Semantic CCX lowering runs before generic decomposition. With ancillas
        // barred the pass is a no-op, so DecomposePass emits the seven-T circuit.
        if (ccx_options.allow_ancilla &&
            (ccx_options.strategy == NWQEC::CCXStrategy::Jones4T ||
             ccx_options.strategy == NWQEC::CCXStrategy::SharedAND)) {
            passes.insert(passes.begin(), NWQEC::PassType::CCX_LOWERING);
        }

        // MCX lowering runs first: the terminal Toffolis it emits are then realized by
        // the CCX pass above, so its choice governs their cost.
        if (mcx_options.strategy != NWQEC::MCXStrategy::Preserve) {
            passes.insert(passes.begin(), NWQEC::PassType::MCX_LOWERING);
        }

        // Add cleanup passes if requested
        if (remove_pauli) {
            passes.push_back(NWQEC::PassType::REMOVE_PAULI);
        }
        
        return transpiler.execute_passes(std::move(circuit_copy), passes, config);
    }
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "NWQEC Python bindings";

#ifdef NWQEC_WITH_GRIDSYNTH_CPP
    m.attr("WITH_GRIDSYNTH_CPP") = py::bool_(NWQEC_WITH_GRIDSYNTH_CPP);
#else
    m.attr("WITH_GRIDSYNTH_CPP") = py::bool_(false);
#endif

    // Proxy returned by Circuit.c_if(...). Gate methods append the gate under the
    // recorded condition and return the proxy, so one c_if chain maps to exactly one
    // emitted `if (...) { ... }` block.
    py::class_<ConditionalBuilder>(m, "ConditionalBuilder")
        .def("__repr__", [](const ConditionalBuilder &b)
             {
                 std::string r = "<ConditionalBuilder bits=[";
                 for (size_t i = 0; i < b.condition.bits.size(); ++i)
                     r += (i ? "," : "") + std::to_string(b.condition.bits[i]);
                 return r + "] value=" + std::to_string(b.condition.value) + ">"; })
#define NWQEC_CB_1Q(NAME, TYPE)                                                       \
    .def(NAME, [](ConditionalBuilder &b, size_t q) -> ConditionalBuilder &            \
         { return b.guarded([&](NWQEC::Circuit &c)                                    \
                            { c.add_operation({NWQEC::Operation::TYPE, {q}}); }); },  \
         py::arg("q"))
#define NWQEC_CB_2Q(NAME, TYPE)                                                            \
    .def(NAME, [](ConditionalBuilder &b, size_t q0, size_t q1) -> ConditionalBuilder &     \
         { return b.guarded([&](NWQEC::Circuit &c)                                         \
                            { c.add_operation({NWQEC::Operation::TYPE, {q0, q1}}); }); },  \
         py::arg("q0"), py::arg("q1"))
#define NWQEC_CB_ROT(NAME, TYPE)                                                              \
    .def(NAME, [](ConditionalBuilder &b, size_t q, double th) -> ConditionalBuilder &         \
         { return b.guarded([&](NWQEC::Circuit &c)                                            \
                            { c.add_operation({NWQEC::Operation::TYPE, {q}, {th}}); }); },    \
         py::arg("q"), py::arg("theta"))
            NWQEC_CB_1Q("x", Type::X) NWQEC_CB_1Q("y", Type::Y) NWQEC_CB_1Q("z", Type::Z)
            NWQEC_CB_1Q("h", Type::H) NWQEC_CB_1Q("s", Type::S) NWQEC_CB_1Q("sdg", Type::SDG)
            NWQEC_CB_1Q("t", Type::T) NWQEC_CB_1Q("tdg", Type::TDG)
            NWQEC_CB_1Q("sx", Type::SX) NWQEC_CB_1Q("sxdg", Type::SXDG)
            NWQEC_CB_2Q("cx", Type::CX) NWQEC_CB_2Q("cz", Type::CZ) NWQEC_CB_2Q("swap", Type::SWAP)
            NWQEC_CB_ROT("rx", Type::RX) NWQEC_CB_ROT("ry", Type::RY) NWQEC_CB_ROT("rz", Type::RZ)
#undef NWQEC_CB_1Q
#undef NWQEC_CB_2Q
#undef NWQEC_CB_ROT
        .def("ccx", [](ConditionalBuilder &b, size_t q0, size_t q1, size_t q2) -> ConditionalBuilder &
             { return b.guarded([&](NWQEC::Circuit &c)
                                { c.add_operation({NWQEC::Operation::Type::CCX, {q0, q1, q2}}); }); },
             py::arg("q0"), py::arg("q1"), py::arg("q2"))
        .def("mcx", [](ConditionalBuilder &b, const std::vector<size_t> &controls, size_t target) -> ConditionalBuilder &
             { return b.guarded([&](NWQEC::Circuit &c)
                                {
                                    std::vector<size_t> qs(controls);
                                    qs.push_back(target);
                                    c.add_operation({NWQEC::Operation::Type::MCX, std::move(qs)}); }); },
             py::arg("controls"), py::arg("target"))
        .def("circuit", [](ConditionalBuilder &b) -> NWQEC::Circuit &
             { return *b.circuit; }, py::return_value_policy::reference,
             "Return the underlying circuit, ending the conditional chain.");

    // Circuit class (owned by Python via unique_ptr)
    py::class_<NWQEC::Circuit, std::unique_ptr<NWQEC::Circuit>>(m, "Circuit")
        // Circuit constructor
        .def(py::init([](size_t num_qubits)
                      {
                 auto c = std::make_unique<NWQEC::Circuit>();
                 if (num_qubits > 0) c->add_qreg("q", num_qubits);
                 return c; }),
             py::arg("num_qubits"))
        .def("x", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             {
                 if (circuit_has_pauli_ops(c))
                     throw std::runtime_error("Cannot mix Pauli-based operations with standard gates in one circuit (PBC-only).");
                 c.add_operation({NWQEC::Operation::Type::X, {q}});
                 return c; }, py::arg("q"), "Apply Pauli-X to qubit q.")
        .def("y", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             {
                 if (circuit_has_pauli_ops(c))
                     throw std::runtime_error("Cannot mix Pauli-based operations with standard gates in one circuit (PBC-only).");
                 c.add_operation({NWQEC::Operation::Type::Y, {q}});
                 return c; }, py::arg("q"), "Apply Pauli-Y to qubit q.")
        .def("z", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             {
                 if (circuit_has_pauli_ops(c))
                     throw std::runtime_error("Cannot mix Pauli-based operations with standard gates in one circuit (PBC-only).");
                 c.add_operation({NWQEC::Operation::Type::Z, {q}});
                 return c; }, py::arg("q"), "Apply Pauli-Z to qubit q.")
        .def("h", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::H, {q}}); return c; }, py::arg("q"), "Apply Hadamard to qubit q.")
        .def("s", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::S, {q}}); return c; }, py::arg("q"), "Apply phase S (π/2 about Z) to qubit q.")
        .def("sdg", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SDG, {q}}); return c; }, py::arg("q"), "Apply S† to qubit q.")
        .def("t", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::T, {q}}); return c; }, py::arg("q"), "Apply T (π/4 about Z) to qubit q.")
        .def("tdg", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::TDG, {q}}); return c; }, py::arg("q"), "Apply T† to qubit q.")
        .def("sx", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SX, {q}}); return c; }, py::arg("q"), "Apply √X to qubit q.")
        .def("sxdg", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SXDG, {q}}); return c; }, py::arg("q"), "Apply (√X)† to qubit q.")
        .def("cx", [](NWQEC::Circuit &c, size_t q0, size_t q1) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::CX, {q0, q1}}); return c; }, py::arg("q0"), py::arg("q1"), "Apply CX(control=q0, target=q1).")
        .def("ccx", [](NWQEC::Circuit &c, size_t q0, size_t q1, size_t q2) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::CCX, {q0, q1, q2}}); return c; }, py::arg("q0"), py::arg("q1"), py::arg("q2"), "Apply CCX(control=q0,q1; target=q2).")
        .def("mcx", [](NWQEC::Circuit &c, const std::vector<size_t> &controls, size_t target) -> NWQEC::Circuit &
             {
                 if (controls.empty())
                     throw std::invalid_argument("mcx needs at least one control");
                 std::vector<size_t> qs(controls);
                 qs.push_back(target);
                 c.add_operation({NWQEC::Operation::Type::MCX, std::move(qs)});
                 return c; },
             py::arg("controls"), py::arg("target"),
             "Apply a multi-controlled X. Lowered by the MCX pass; see mcx_strategy.")
        .def("cz", [](NWQEC::Circuit &c, size_t q0, size_t q1) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::CZ, {q0, q1}}); return c; }, py::arg("q0"), py::arg("q1"), "Apply CZ between q0 and q1.")
        .def("swap", [](NWQEC::Circuit &c, size_t q0, size_t q1) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SWAP, {q0, q1}}); return c; }, py::arg("q0"), py::arg("q1"), "Swap states of q0 and q1.")
        .def("rx", [](NWQEC::Circuit &c, size_t q, double th) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RX, {q}, {th}}); return c; }, py::arg("q"), py::arg("theta"))
        .def("rxp", [](NWQEC::Circuit &c, size_t q, double x) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RX, {q}, {x * M_PI}}); return c; }, py::arg("q"), py::arg("x_pi"))
        .def("ry", [](NWQEC::Circuit &c, size_t q, double th) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RY, {q}, {th}}); return c; }, py::arg("q"), py::arg("theta"))
        .def("ryp", [](NWQEC::Circuit &c, size_t q, double x) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RY, {q}, {x * M_PI}}); return c; }, py::arg("q"), py::arg("x_pi"))
        .def("rz", [](NWQEC::Circuit &c, size_t q, double th) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RZ, {q}, {th}}); return c; }, py::arg("q"), py::arg("theta"))
        .def("rzp", [](NWQEC::Circuit &c, size_t q, double x) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RZ, {q}, {x * M_PI}}); return c; }, py::arg("q"), py::arg("x_pi"))
        .def("measure", [](NWQEC::Circuit &c, size_t q, size_t b) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::MEASURE, {q}, {}, {b}}); return c; }, py::arg("q"), py::arg("cbit"))
        .def("reset", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RESET, {q}}); return c; }, py::arg("q"))
        .def("barrier", [](NWQEC::Circuit &c, const std::vector<size_t> &qs) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::BARRIER, qs}); return c; }, py::arg("qubits"))
        // Clean Pauli helpers: accept only a string
        .def("t_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::T_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Apply rotation by π/4 about the given Pauli string (e.g., '+XIZ').")
        .def("m_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::M_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Measure the given multi‑qubit Pauli string (projective measurement).")
        .def("s_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::S_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Apply rotation by π/2 about the given Pauli string.")
        .def("z_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::Z_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Apply rotation by π about the given Pauli string.")
        .def("c_if", [](NWQEC::Circuit &c, const std::vector<size_t> &cbits, uint64_t value)
             {
                 if (cbits.empty())
                     throw std::runtime_error("c_if requires at least one classical bit.");
                 return ConditionalBuilder{&c, NWQEC::ClassicalCondition(cbits, value)}; },
             py::arg("cbits"), py::arg("value") = 1,
             "Begin a classically conditioned block. Gates chained onto the returned "
             "builder execute only when the given classical bits equal `value`.\n"
             "Example: circuit.c_if([0], 1).cz(0, 1)")
        .def("num_qubits", &NWQEC::Circuit::get_num_qubits)
        .def("num_bits", &NWQEC::Circuit::get_num_bits)
        .def("has_measurement", &NWQEC::Circuit::has_measurement)
        .def("has_feedforward", &NWQEC::Circuit::has_feedforward)
        .def("is_dynamic", &NWQEC::Circuit::is_dynamic)
    .def("count_ops", &NWQEC::Circuit::count_ops)
    .def("is_clifford_t", &NWQEC::Circuit::is_clifford_t)
        .def("stats", &circuit_stats)
        .def("duration", &NWQEC::Circuit::duration, py::arg("code_distance"))
        .def("depth", &NWQEC::Circuit::depth)
        .def("to_qasm", &circuit_to_qasm)
        .def("to_qasm_str", &circuit_to_qasm)
        .def("save_qasm", &circuit_save_qasm, py::arg("path"))
        .def("to_qasm_file", &circuit_save_qasm, py::arg("filename"));

    // Module-level transforms: clean entrypoints
    m.def(
        "to_clifford_t",
        [](const NWQEC::Circuit &circuit, bool keep_ccx, const std::string &rz_err, py::object epsilon,
           const std::string &ccx_strategy, size_t max_ancillas, bool allow_ancilla,
           const std::string &mcx_strategy, size_t mcx_max_ancillas)
        {
            auto ccx_options = resolve_ccx_options(ccx_strategy, keep_ccx, max_ancillas, allow_ancilla);
            auto mcx_options = resolve_mcx_options(mcx_strategy, mcx_max_ancillas, allow_ancilla);
            return apply_transforms(circuit,
                                    /*to_pbc=*/false,
                                    /*to_clifford_reduction=*/false,
                                    /*keep_cx=*/false,
                                    /*t_pauli_opt=*/false,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/keep_ccx,
                                    /*silent=*/true,
                                    rz_err,
                                    epsilon,
                                    ccx_options,
                                    mcx_options);
        },
        py::arg("circuit"),
        py::arg("keep_ccx") = false,
        py::arg("rz_err") = "per-gate",
        py::arg("epsilon") = py::none(),
        py::arg("ccx_strategy") = "",
        py::arg("max_ancillas") = 0,
        py::arg("allow_ancilla") = true,
        py::arg("mcx_strategy") = "",
        py::arg("mcx_max_ancillas") = 0,
        "Convert the input circuit to a Clifford+T-only circuit and return a new Circuit.\n"
        "- keep_ccx: preserve CCX gates during decomposition\n"
        "- rz_err: RZ synthesis error policy: 'per-gate', 'total', or 'relative'\n"
        "- epsilon: optional value for the selected RZ error policy\n"
        "- ccx_strategy: 'standard-7t' (7 T per CCX, no ancilla), 'jones-4t' (4 T per CCX,\n"
        "  one clean ancilla plus a mid-circuit measurement and a conditioned CZ),\n"
        "  'shared-and' (4 T per GROUP of CCX gates sharing a control pair), or\n"
        "  'preserve' to leave CCX gates intact\n"
        "- max_ancillas: ancilla budget for 'shared-and'; 0 means unlimited, i.e. as many\n"
        "  as the grouping requires. Not used by 'jones-4t', whose groups are always of\n"
        "  size one and therefore never overlap\n"
        "- allow_ancilla: when False no ancilla may be used, so the ancilla-free 7 T\n"
        "  decomposition is emitted whatever ccx_strategy names. This is a property of\n"
        "  the target rather than a preference, so it overrides the strategy\n"
        "- mcx_strategy: how each MCX is realized, all with clean ancillae: 'and-ladder'\n"
        "  and 'and-tree' (n-1 ancillas, 4n-4 T), 'vchain-m15' (n-2 ancillas, 4n-1 T),\n"
        "  'hybrid-budget' (A ancillas, 8n-4A-5 T), 'kg24-1-clean' (1 ancilla, 8n-9 T),\n"
        "  'kg24-2-clean' (2 ancillas, 8n-9 T, log depth), or 'preserve' (default)\n"
        "- mcx_max_ancillas: MCX ancilla budget, and the budget A for 'hybrid-budget';\n"
        "  0 means unlimited, which for 'hybrid-budget' means A = n-1");

    m.def(
        "to_pbc",
        [](const NWQEC::Circuit &circuit, bool keep_cx, bool optimize_t_count, const std::string &rz_err, py::object epsilon,
           const std::string &ccx_strategy, size_t max_ancillas, bool allow_ancilla,
           const std::string &mcx_strategy, size_t mcx_max_ancillas)
        {
            auto ccx_options = resolve_ccx_options(ccx_strategy, /*keep_ccx=*/false, max_ancillas, allow_ancilla);
            auto mcx_options = resolve_mcx_options(mcx_strategy, mcx_max_ancillas, allow_ancilla);
            // Use optimized PBC pipeline if T-optimization is requested
            if (optimize_t_count) {
                NWQEC::Transpiler transpiler;
                NWQEC::PassConfig config;
                config.keep_cx = keep_cx;
                config.silent = true;
                config.ccx_lowering = ccx_options;
                config.mcx_lowering = mcx_options;
                configure_rz_error(config, rz_err, epsilon);
                
                auto circuit_copy = std::make_unique<NWQEC::Circuit>(circuit);
                auto passes = NWQEC::PassSequences::TO_PBC_OPTIMIZED;
                if (ccx_options.allow_ancilla &&
                    (ccx_options.strategy == NWQEC::CCXStrategy::Jones4T ||
                     ccx_options.strategy == NWQEC::CCXStrategy::SharedAND)) {
                    passes.insert(passes.begin(), NWQEC::PassType::CCX_LOWERING);
                }
                if (mcx_options.strategy != NWQEC::MCXStrategy::Preserve) {
                    passes.insert(passes.begin(), NWQEC::PassType::MCX_LOWERING);
                }

                return transpiler.execute_passes(std::move(circuit_copy), passes, config);
            } else {
                return apply_transforms(circuit,
                                        /*to_pbc=*/true,
                                        /*to_clifford_reduction=*/false,
                                        /*keep_cx=*/keep_cx,
                                        /*t_pauli_opt=*/false,
                                        /*remove_pauli=*/false,
                                        /*keep_ccx=*/false,
                                        /*silent=*/true,
                                        rz_err,
                                        epsilon,
                                        ccx_options,
                                        mcx_options);
            }
        },
        py::arg("circuit"),
        py::arg("keep_cx") = false,
        py::arg("optimize_t_count") = false,
        py::arg("rz_err") = "per-gate",
        py::arg("epsilon") = py::none(),
        py::arg("ccx_strategy") = "",
        py::arg("max_ancillas") = 0,
        py::arg("allow_ancilla") = true,
        py::arg("mcx_strategy") = "",
        py::arg("mcx_max_ancillas") = 0,
        "Transpile the input circuit to a Pauli-Based Circuit (PBC) form and return a new Circuit.\n"
        "- keep_cx: preserve CX gates where possible in the PBC form\n"
        "- optimize_t_count: apply T-count optimization after PBC conversion\n"
        "- rz_err: RZ synthesis error policy: 'per-gate', 'total', or 'relative'\n"
        "- epsilon: optional value for the selected RZ error policy");

    m.def(
        "to_clifford_reduction",
        [](const NWQEC::Circuit &circuit, const std::string &rz_err, py::object epsilon)
        {
            return apply_transforms(circuit,
                                    /*to_pbc=*/false,
                                    /*to_clifford_reduction=*/true,
                                    /*keep_cx=*/false,
                                    /*t_pauli_opt=*/false,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/false,
                                    /*silent=*/true,
                                    rz_err,
                                    epsilon);
        },
        py::arg("circuit"),
        py::arg("rz_err") = "per-gate",
        py::arg("epsilon") = py::none(),
        "Apply the Clifford reduction optimization and return a new Circuit.\n"
        "This optimization preserves circuit parallelism while reducing non-T overhead.\n"
        "Based on the technique from: Wang et al. 'Optimizing FTQC Programs through QEC Transpiler and Architecture Codesign' (2024)\n"
        "- rz_err: RZ synthesis error policy: 'per-gate', 'total', or 'relative'\n"
        "- epsilon: optional value for the selected RZ error policy");


    // fuse_t: apply only the T-Pauli fusion stage within the PBC pipeline
    m.def(
        "fuse_t",
        [](const NWQEC::Circuit &circuit, const std::string &rz_err, py::object epsilon)
        {
            return apply_transforms(circuit,
                                    /*to_pbc=*/false,
                                    /*to_clifford_reduction=*/false,
                                    /*keep_cx=*/false,
                                    /*t_pauli_opt=*/true,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/false,
                                    /*silent=*/true,
                                    rz_err,
                                    epsilon);
        },
        py::arg("circuit"),
        py::arg("rz_err") = "per-gate",
        py::arg("epsilon") = py::none(),
        "Optimize the number of T rotations within a Pauli-Based Circuit (PBC) and return a new Circuit.\n"
        "- rz_err: RZ synthesis error policy: 'per-gate', 'total', or 'relative'\n"
        "- epsilon: optional value for any RZ synthesis still required");

    m.def(
        "get_clifford_t_counts",
        [](const NWQEC::Circuit &circuit, bool keep_ccx, const std::string &rz_err, py::object epsilon)
        {
            NWQEC::RzErrorPolicy policy = parse_rz_err_policy(rz_err);
            double resolved_epsilon = epsilon.is_none()
                                          ? NWQEC::default_rz_error_epsilon(policy)
                                          : epsilon.cast<double>();

            return NWQEC::get_clifford_t_counts(circuit, policy, resolved_epsilon, keep_ccx);
        },
        py::arg("circuit"),
        py::arg("keep_ccx") = false,
        py::arg("rz_err") = "per-gate",
        py::arg("epsilon") = py::none(),
        "Return exact Clifford+T gate counts without generating the final circuit.\n"
        "- keep_ccx: preserve CCX gates during decomposition\n"
        "- rz_err: RZ synthesis error policy: 'per-gate', 'total', or 'relative'\n"
        "- epsilon: optional value for the selected RZ error policy");

    m.def(
        "mcx_cost",
        [](size_t num_controls, const std::string &strategy, size_t budget)
        {
            auto cost = NWQEC::MCXLoweringPass::mcx_cost(num_controls,
                                                         parse_mcx_strategy(strategy), budget);
            py::dict d;
            d["and_count"] = cost.and_count;
            d["rccx_count"] = cost.rccx_count;
            d["ccx_count"] = cost.ccx_count;
            d["clean_ancillas"] = cost.clean_ancillas;
            d["result_bits"] = cost.result_bits;
            d["depth_class"] = std::string(cost.depth_class);
            d["t_standard7t"] = cost.t_standard7t();
            d["t_jones4t"] = cost.t_jones4t();
            d["ancillas_jones4t"] = cost.ancillas_jones4t();
            return d;
        },
        py::arg("num_controls"), py::arg("strategy"), py::arg("budget") = 0,
        "Static resource model for one C^nX under the named MCX strategy.\n"
        "Returns gate-class counts plus the T-count under each downstream CCX strategy,\n"
        "since the terminal Toffolis are realized by the CCX pass rather than here.");

    m.def("load_qasm", [](const std::string &filename)
          {
        NWQEC::QASMParser p;
        if (!p.parse_file(filename))
        {
            throw std::runtime_error("Failed to parse QASM: " + p.get_error_message());
        }
        return p.get_circuit(); }, py::arg("filename"));
}
