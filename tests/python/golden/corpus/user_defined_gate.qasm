OPENQASM 2.0;
include "qelib1.inc";
gate mygate a,b {
  h a;
  cx a,b;
  t b;
}
gate wrapper a,b,c {
  mygate a,b;
  ccx a,b,c;
}
qreg q[3];
creg c[3];
h q[0];
mygate q[0],q[1];
wrapper q[0],q[1],q[2];
