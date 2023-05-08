function [omega_p] = A_matrix(quantities, omega)

  T1 = quantities.T1; U1 = quantities.U1; T2 = quantities.T2; U2 = quantities.U2;
  F11 = quantities.F11; F12 = quantities.F12; F21 = quantities.F21; F22 = quantities.F22;
  G11 = quantities.G11; G12 = quantities.G12; G21 = quantities.G21; G22 = quantities.G22;
  H11 = quantities.H11; H21 = quantities.H21;

  syms A

  syms omega_p_sym
  A(1,1) = -(omega_p_sym^2)*T1 + U1 - i*omega(1)*G11 - F11 - omega(1)*H11;
  A(1,2) = -F12 - i*omega(1)*G12;
  A(2,1) = -F21 - i*omega(1)*G21 - omega(1)*H21;
  A(2,2) = -(omega(1)^2)*T2 + U2 - i*omega(1)*G22 - F22;

  charpoly = det(A) == 0;
  B = double(solve(charpoly));
  omega_p(1) = B(real(B)>0);

  clear omega_p_sym; syms omega_p_sym
  A(1,1) = -(omega(2)^2)*T1 + U1 - i*omega(2)*G11 - F11 - omega(2)*H11;
  A(1,2) = -F12 - i*omega(2)*G12;
  A(2,1) = -F21 - i*omega(2)*G21 - omega(2)*H21;
  A(2,2) = -(omega_p_sym^2)*T2 + U2 - i*omega(2)*G22 - F22;

  charpoly = det(A) == 0;
  B = double(solve(charpoly));
  omega_p(2) = B(real(B)>0);

end
