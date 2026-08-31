OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];

t_pauli +ZZI;
t_pauli -YYI;
t_pauli +YYX;
t_pauli +ZZX;
t_pauli -YYI;
t_pauli -XXX;
t_pauli +XXI;
t_pauli -ZZI;
t_pauli +IIX;
m_pauli +XXI;
m_pauli -YYI;
m_pauli +IIZ;
