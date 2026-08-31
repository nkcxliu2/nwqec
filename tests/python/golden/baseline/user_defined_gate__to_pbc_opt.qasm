OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

t_pauli +YYX;
t_pauli +ZZX;
t_pauli -XXX;
t_pauli +XXI;
t_pauli +IIX;
m_pauli +IIZ;
m_pauli -YYI;
m_pauli +XXI;
