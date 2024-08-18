function [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,all_log_eigs_stiff_mat)
% Given a minimizer "zeq" with "nbp" base pairs, and precomputed array
%   of log of eigenvalues of stiffness matrix, finish the Laplace 
%   approximation by computing energy "eneq", log-base-10 of det in 
%   formula, and contribution to j in nM
  nframes = nbp+1;
  [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess_nopen(zeq);
  hess_wrt_k = convert_to_hess_wrt_k(zeq,hesseq);
  newzlen=6*nframes+6*(nbp-1);
  small_eigs=eigs(hess_wrt_k,newzlen/2,'SM');
  large_eigs=eigs(hess_wrt_k,newzlen/2,'LM');
  large_eigs=large_eigs(newzlen/2:-1:1);   % Sort largest eigs in incr order
  all_log_eigs_hess_wrt_k=[log10(small_eigs);log10(large_eigs)];

% Compute log of product of gammas
  sum_of_log_gamma = compute_sum_log_gamma(zeq);
  log10_overall_det = sum(all_log_eigs_stiff_mat)/2-sum(all_log_eigs_hess_wrt_k)/2;
  log10_overall_det = log10_overall_det+3+sum_of_log_gamma*2;
  %term1_is = sum_of_log_gamma*2
  %term2_is = sum(all_log_eigs_stiff_mat)/2
  %term3_is = sum(all_log_eigs_hess_wrt_k)/2
  jfac=10^(log10_overall_det)*exp(-eneq)/(2*pi)^3*pi^2;  % Length in angstroms
  jfac_in_mol = jfac*10^4/6;
  jfac_cont = jfac_in_mol*10^9;

  return
end