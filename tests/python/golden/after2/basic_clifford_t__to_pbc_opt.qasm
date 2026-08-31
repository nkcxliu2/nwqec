OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

t_pauli +IZI;
t_pauli -IIZ;
m_pauli -IIY;
m_pauli -XZI;
m_pauli +XII;
