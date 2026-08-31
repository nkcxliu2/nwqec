OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

t_pauli +XII;
m_pauli +IIZ;
m_pauli +XZI;
m_pauli +XII;
