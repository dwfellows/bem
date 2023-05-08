function [omega_soln] = omega_iterative(quantities, omega0)

  T1 = quantities.T1; U1 = quantities.U1; T2 = quantities.T2; U2 = quantities.U2;
  omega0 = [sqrt(U1/T1), sqrt(U2/T2)];

  chi = 0.5; % Relaxation coefficient
  eps1 = 10^(-3);
  eps2 = 10^(-3);

  converged = 0;
  while (converged < 1)

    omega_int = A_matrix(quantities, omega0);
    omega_new = (1-chi)*omega_int + chi*omega0;

    check = abs(((omega_new - omega0)./omega0));
    if (check < eps1)
      converged = 1;
      omega_soln = omega_new;
    else
      omega0 = omega_new;
    end

  end

end
