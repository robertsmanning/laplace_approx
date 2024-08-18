function jfac_in_nmol = approx_j_multi(nbpmin,nbpmax)
% Compute Laplace approx to J for rigid base pair model for a range of DNA
% lengths

% Inputs: nbpmin = smallest length (in bp), nbpmax = largest length (in bp)
%    (DNA model params set on lines 16-25 below

% Outputs: jfac_in_nmol = J in nM for last DNA considered
%    (Results for all mols stored in output file called "results_multi")

fid_out = fopen("results_multi","w");
for nbp = nbpmin:nbpmax
  % Add phantom frame
  nframes = nbp+1;
  % Shape and stiffness params, including seq-dependent theta1hat, theta2hat
  global stiffmat
  B1 = 1/(3.4/570); B3 = 2.0*B1/3.0; A1 = B1; A3 = 10*A1;
  stiffmat = diag([B1 B1 B3 A1 A1 A3]);

  global whats
  theta3hat = 2*pi/10.5; zeta3hat = 3.4;
  whats = zeros(6,nbp); whats(3,:) = theta3hat; whats(6,:) = zeta3hat;
  for i = 18:122
    whats(1,i) = (pi/2)/105*cos(theta3hat*(i-18));
    whats(2,i) = (pi/2)/105*sin(-theta3hat*(i-18));
  end

  % Value of q4 (+1 or -1) at far end of DNA
  global q4_at_1

  % Compute sum of log of eigs of stiffness matrix
  sum_of_log_eigs_stiff = nbp*(2*log10(A1)+log10(A3)+2*log10(B1)+log10(B3));

  dt = 1/nbp; R = (zeta3hat*nbp)/(2*pi);
  lk0 = round(nbp*theta3hat/(2*pi)); 

  jfac_in_nmol = 0.0; summary_mat = [];

  %  Find contribution from Lk = Lk0 (use best register to build init guess)
  w = 2*pi*lk0; lk = lk0;
  enmin = 9999; tht = 0.0;
  for ireg = 0:11
    enreg = energy_of_twisted_circle(nbp,ireg*pi/6,lk);
    if enreg < enmin
      enmin = enreg; tht = ireg*pi/6;
    end
  end
  [rs,qs,zvec] = build_twisted_circle(nbp,tht,w);
  q4_at_1 = qs(4,nframes);
  fprintf("\nSearching for energy minimizer, %d bp, %d init Lk\n",nbp,lk);
  [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
  lk_comp = compute_lk(nbp,zeq);
  max_normq_err = compute_max_normq_err(nbp,zeq);
  fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

  [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
  fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));

  [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,sum_of_log_eigs_stiff);
  jfac_in_nmol = jfac_in_nmol + jfac_cont;
  summary_mat = [summary_mat; lk_comp log10(exp(-eneq)) log10_overall_det-log10(8*pi*6)+13 all_time];
  
  % Consider Lk = Lk0-1
  w = 2*pi*(lk0-1); lk = lk0-1;
  enmin = 9999; tht = 0.0;
  for ireg = 0:11
    enreg = energy_of_twisted_circle(nbp,ireg*pi/6,lk);
    if enreg < enmin
      enmin = enreg; tht = ireg*pi/6;
    end
  end
  [rs,qs,zvec] = build_twisted_circle(nbp,tht,w);
  q4_at_1 = qs(4,nframes);
  fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
  [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
  lk_comp = compute_lk(nbp,zeq);
  max_normq_err = compute_max_normq_err(nbp,zeq);
  fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

  [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
  fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));

  [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,sum_of_log_eigs_stiff);
  jfac_in_nmol = jfac_in_nmol + jfac_cont;
  summary_mat = [summary_mat; lk_comp log10(exp(-eneq)) log10_overall_det-log10(8*pi*6)+13 all_time];
  
  % Consider Lk = Lk0-2 only if Lk0-1 had no strand passage and decreased E
  if abs(lk_comp-(lk0-1)) > 0.5 || summary_mat(2,2)<summary_mat(1,2)
    fprintf("\nNot searching for energy minimizer with %d init Lk\n",lk0-2);
    summary_mat = [summary_mat; lk0-2 -99.0 99.0 0.0];
  else
    w = 2*pi*(lk0-2); lk = lk0-2;
    enmin = 9999; tht = 0.0;
    for ireg = 0:11
      enreg = energy_of_twisted_circle(nbp,ireg*pi/6,lk);
      if enreg < enmin
        enmin = enreg; tht = ireg*pi/6;
      end
    end
    [rs,qs,zvec] = build_twisted_circle(nbp,tht,w);
    q4_at_1 = qs(4,nframes);
    fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
    [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
    lk_comp = compute_lk(nbp,zeq);
    max_normq_err = compute_max_normq_err(nbp,zeq);
    fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

    [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
    fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));

    [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,sum_of_log_eigs_stiff);
    
    is_new = 1;
    for i2 = 1:2
      if abs(lk_comp-summary_mat(i2,1)) < 0.5 && abs(log10(exp(-eneq))-summary_mat(i2,2)) < 0.1
        is_new = 0;
      end
    end
    if is_new == 1
      jfac_in_nmol = jfac_in_nmol + jfac_cont;
      summary_mat = [summary_mat; lk_comp log10(exp(-eneq)) log10_overall_det-log10(8*pi*6)+13 all_time];
    else
      summary_mat = [summary_mat; lk_comp -99.0 99.0 all_time];
    end
  end

  % Consider Lk = Lk0+1 unless E(Lk0) > E(Lk0-1) or link of E(Lk0-1) is Lk0+1
  if summary_mat(1,2) < summary_mat(2,2) && abs(summary_mat(2,1)-(lk0-1)) < 0.5
    summary_mat = [summary_mat; lk0+1 -99.0 99.0 0.0];
    fprintf("\nNot searching with %d init Lk since E would be too high\n",lk0+1);
  else
    if abs(summary_mat(2,1)-(lk0+1)) < 0.5
      summary_mat = [summary_mat; lk0+1 -99.0 99.0 0.0];
      fprintf("\nNot searching with %d init Lk since already found\n",lk0+1);
    else
      w = 2*pi*(lk0+1); lk = lk0+1;
      enmin = 9999; tht = 0.0;
      for ireg = 0:11
        enreg = energy_of_twisted_circle(nbp,ireg*pi/6,lk);
        if enreg < enmin
          enmin = enreg; tht = ireg*pi/6;
        end
      end
      [rs,qs,zvec] = build_twisted_circle(nbp,tht,w);
      q4_at_1 = qs(4,nframes);
      fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
      [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
      lk_comp = compute_lk(nbp,zeq);
      max_normq_err = compute_max_normq_err(nbp,zeq);
      fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

      [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
      fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));

      [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,sum_of_log_eigs_stiff);
    
      is_new = 1;
      for i2 = 1:3
        if abs(lk_comp-summary_mat(i2,1)) < 0.5 && abs(log10(exp(-eneq))-summary_mat(i2,2)) < 0.1
          is_new = 0;
        end
      end
      if is_new == 1
        jfac_in_nmol = jfac_in_nmol + jfac_cont;
        summary_mat = [summary_mat; lk_comp log10(exp(-eneq)) log10_overall_det-log10(8*pi*6)+13 all_time];
      else
        summary_mat = [summary_mat; lk_comp -99.0 99.0 all_time];
      end
    end
  end

  % Consider Lk = Lk0+2 only if Lk0+1 reduced energy
  if (summary_mat(4,2) > summary_mat(1,2) && abs(summary_mat(4,1)-(lk0+1)) < 0.5) || (summary_mat(2,2) > summary_mat(1,2) && abs(summary_mat(2,1)-(lk0+1))<0.5)
    w = 2*pi*(lk0+2); lk = lk0+2;
    enmin = 9999; tht = 0.0;
    for ireg = 0:11
      enreg = energy_of_twisted_circle(nbp,ireg*pi/6,lk);
      if enreg < enmin
        enmin = enreg; tht = ireg*pi/6;
      end
    end
    [rs,qs,zvec] = build_twisted_circle(nbp,tht,w);
    q4_at_1 = qs(4,nframes);
    fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
    [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
    lk_comp = compute_lk(nbp,zeq);
    max_normq_err = compute_max_normq_err(nbp,zeq);
    fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

    [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
    fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));

    [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,sum_of_log_eigs_stiff);
    
    is_new = 1;
    for i2 = 1:4
      if abs(lk_comp-summary_mat(i2,1)) < 0.5 && abs(log10(exp(-eneq))-summary_mat(i2,2)) < 0.1
        is_new = 0;
      end
    end
    if is_new == 1
      jfac_in_nmol = jfac_in_nmol + jfac_cont;
      summary_mat = [summary_mat; lk_comp log10(exp(-eneq)) log10_overall_det-log10(8*pi*6)+13 all_time];
    else
      summary_mat = [summary_mat; lk_comp -99.0 99.0 all_time];
    end
  else
    summary_mat = [summary_mat; lk0+2 -99.0 99.0 0.0];
    fprintf("\nNot searching for energy minimizer with %d init Lk\n",lk0+2);
  end

  summary_mat = sortrows(summary_mat,2,'descend');
  log_jfac_in_nmol = log10(jfac_in_nmol);
  fprintf("Overall, J = %6.3f, log10J = %6.3f\n",jfac_in_nmol,log_jfac_in_nmol);
  fprintf(fid_out,"%d %8.5f ",nbp,log_jfac_in_nmol);
  for i = 1:5
      fprintf(fid_out,"%d %6.3f %6.3f %7.1f   ",round(summary_mat(i,1)),summary_mat(i,2),summary_mat(i,3),summary_mat(i,4));
  end
  fprintf(fid_out,"\n");
end
fclose(fid_out)
