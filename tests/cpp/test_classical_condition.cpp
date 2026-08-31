// Tests for classically conditioned operations that are not reachable from Python:
// DAG classical dependency edges, and IR-level condition handling.
#include "nwqec/core/circuit.hpp"
#include "nwqec/core/dag_circuit.hpp"
#include "nwqec/core/operation.hpp"

#include <iostream>
#include <string>

static int failures = 0;

static void check(const std::string &name, bool ok)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name << "\n";
    if (!ok)
        ++failures;
}

static void test_condition_ir()
{
    std::cout << "[IR]\n";
    NWQEC::Operation op(NWQEC::Operation::Type::CZ, {0, 1});
    check("gate is unconditional by default", !op.is_conditional());

    NWQEC::Operation guarded = op.with_condition(NWQEC::ClassicalCondition({3}, 1));
    check("with_condition does not mutate the original", !op.is_conditional());
    check("with_condition returns a conditional copy", guarded.is_conditional());
    check("condition bits survive", guarded.get_condition()->bits == std::vector<size_t>{3});

    NWQEC::Operation copied = guarded;
    check("condition survives copy construction", copied.is_conditional());
    NWQEC::Operation assigned(NWQEC::Operation::Type::X, {0});
    assigned = guarded;
    check("condition survives copy assignment", assigned.is_conditional());

    bool threw = false;
    try
    {
        NWQEC::Operation m(NWQEC::Operation::Type::MEASURE, {0}, {}, {0});
        m.set_condition(NWQEC::ClassicalCondition({1}, 1));
    }
    catch (const std::exception &)
    {
        threw = true;
    }
    check("conditioning a MEASURE is rejected", threw);

    threw = false;
    try
    {
        NWQEC::ClassicalCondition bad({0}, 5); // 5 does not fit in one bit
        (void)bad;
    }
    catch (const std::exception &)
    {
        threw = true;
    }
    check("value wider than the bit vector is rejected", threw);
}

static void test_dag_classical_edges()
{
    std::cout << "[DAG]\n";
    NWQEC::DAGCircuit dag;
    dag.add_qreg("q", 3);
    dag.add_creg("c", 1);

    dag.add_operation(NWQEC::Operation(NWQEC::Operation::Type::H, {0}));                // idx 0
    dag.add_operation(NWQEC::Operation(NWQEC::Operation::Type::MEASURE, {0}, {}, {0})); // idx 1
    NWQEC::Operation cz(NWQEC::Operation::Type::CZ, {1, 2});
    cz.set_condition(NWQEC::ClassicalCondition({0}, 1));
    dag.add_operation(cz); // idx 2

    // The correction shares no qubit with the measurement, so the only thing that can
    // order them is the classical edge.
    bool found = false;
    for (const auto &pred : dag.get_predecessors(2))
    {
        if (pred.is_classical && pred.node == 1 && pred.qubit == 0)
            found = true;
    }
    check("conditional op depends on the measurement that wrote its bit", found);

    bool succ = false;
    for (const auto &s : dag.get_successors(1))
    {
        if (s.is_classical && s.node == 2)
            succ = true;
    }
    check("measurement lists the conditional op as a successor", succ);

    bool spurious = false;
    for (const auto &pred : dag.get_predecessors(2))
    {
        if (!pred.is_classical && pred.node == 1)
            spurious = true;
    }
    check("no spurious quantum edge is created", !spurious);

    bool threw = false;
    try
    {
        NWQEC::Operation bad(NWQEC::Operation::Type::X, {1});
        bad.set_condition(NWQEC::ClassicalCondition({99}, 1));
        dag.add_operation(bad);
    }
    catch (const std::out_of_range &)
    {
        threw = true;
    }
    check("out-of-range condition bit is rejected", threw);
}

static void test_unconditional_dag_unchanged()
{
    std::cout << "[DAG regression]\n";
    NWQEC::DAGCircuit dag;
    dag.add_qreg("q", 2);
    dag.add_operation(NWQEC::Operation(NWQEC::Operation::Type::H, {0}));
    dag.add_operation(NWQEC::Operation(NWQEC::Operation::Type::CX, {0, 1}));
    check("unconditional circuit still builds quantum edges",
          dag.get_predecessors(1).size() == 1 && !dag.get_predecessors(1)[0].is_classical);
    check("unconditional circuit reports no feedforward", !dag.has_feedforward());
}

int main()
{
    test_condition_ir();
    test_dag_classical_edges();
    test_unconditional_dag_unchanged();
    std::cout << (failures ? "FAILED" : "ALL PASSED") << " (" << failures << " failures)\n";
    return failures ? 1 : 0;
}
