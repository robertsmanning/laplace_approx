function [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,sum_of_log_eigs_stiff_mat)
% Given a local minimizer "zeq" with "nbp" base pairs, and a previously
%  computed sum of the log of the eigenvalues of the stiffness matrix,
%  finish the Laplace calculation by finding the energy, the log-base-10
%  of the determinant, and the contribution to J
  nframes = nbp+1;
  [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess_nopen(zeq);
  hess_wrt_k = convert_to_hess_wrt_k(zeq,hesseq);
  newzlen=6*(nbp-1);
  small_eigs=eigs(hess_wrt_k,newzlen/2,'SM');
  large_eigs=eigs(hess_wrt_k,newzlen/2,'LM');
  large_eigs=large_eigs(newzlen/2:-1:1);   % Sort largest eigs in incr order
  all_log_eigs_hess_wrt_k=[log10(small_eigs);log10(large_eigs)];

% Compute log of product of gammas
  sum_of_log_gamma = compute_sum_log_gamma(zeq);
  log10_overall_det = sum_of_log_eigs_stiff_mat/2-sum(all_log_eigs_hess_wrt_k)/2;
  log10_overall_det = log10_overall_det+3*log10(2)+sum_of_log_gamma*2;
  jfac=10^(log10_overall_det)*exp(-eneq)/(2*pi)^3*pi^2;  % Length in angstroms
  jfac_in_mol = jfac*10^4/6;
  jfac_cont = jfac_in_mol*10^9;

  return
end
