function [ret] = plot_soln(col_pts, p_bem_orig, p_bem_soln, mode_ind, mach_ind, M_list)

  ret = 0; afs = 18; fs = 22;
  figure;
  subplot(211);
  plot(col_pts, real(p_bem_orig(mode_ind,:,mach_ind)), 'r', 'LineWidth', 2); grid on; hold on;
  plot(col_pts, real(p_bem_soln(mode_ind,:,mach_ind)), 'b', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  legend('$\omega_0$', '$\omega_f$', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northeast');
  ylabel('Re( $p^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  subplot(212);
  plot(col_pts, imag(p_bem_orig(mode_ind,:,mach_ind)), 'r', 'LineWidth', 2); grid on; hold on;
  plot(col_pts, imag(p_bem_soln(mode_ind,:,mach_ind)), 'b', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  %legend('$\omega_0$', '$\omega_f$', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northeast');
  ylabel('Im( $p^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['Comparison of Initial and Final BEM Solutions, Mode $\psi_{', num2str(mode_ind), '}, M = ', num2str(M_list(mach_ind)), '$']), 'interpreter', 'latex', 'fontsize', fs);


  ret = 1;

end
