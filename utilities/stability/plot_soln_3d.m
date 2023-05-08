function [ret] = plot_soln(col_pts, p_bem_orig, p_bem_soln, mode_ind, mach_ind, M_list)

  ret = 0; afs = 18; fs = 22;

  figure;
  subplot(221);
  scatter(col_pts(:,1), col_pts(:,2), 24, real(p_bem_orig(mode_ind,:,mach_ind))); grid on; hold on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
  ylabel('Re ( $p / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  title('$\omega_0$', 'interpreter', 'latex', 'fontsize', fs);
  subplot(222);
  scatter(col_pts(:,1), col_pts(:,2), 24, real(p_bem_soln(mode_ind,:,mach_ind))); grid on; hold on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
  title('$\omega_f$', 'interpreter', 'latex', 'fontsize', fs);
  subplot(223);
  scatter(col_pts(:,1), col_pts(:,2), 24, imag(p_bem_orig(mode_ind,:,mach_ind))); grid on; hold on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
  ylabel('Im ( $p / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  subplot(224);
  scatter(col_pts(:,1), col_pts(:,2), 24, imag(p_bem_soln(mode_ind,:,mach_ind))); grid on; hold on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
  sgtitle(strcat(['Comparison of Initial and Final BEM Solutions, Mode $W_{', num2str(mode_ind), '}, M = ', num2str(M_list(mach_ind)), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
  ret = 1;

end
