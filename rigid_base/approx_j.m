function log_jfac_in_nmol = approx_j(dnaseq)
% Given a filename "dnaseq" that holds a DNA sequence, finds the Laplace
%  approximation to J using the cgDNA model
%  (lcvmwww.epfl.ch/research/cgDNA/downloads.php#cgDNA)

% Shape and stiffness parameters from cgDNA
global whats
global mmats
global nmats
global omats
global pmats
global qmats
% Value of q4 (+1 or -1) at far end of DNA
global q4_at_1

fid_out = fopen("results_file","w");

  fid = fopen(dnaseq,'r');
  sequence=fscanf(fid,'%s')  
  fclose(fid);
  sequence = strcat(sequence,sequence(1)); % Add the phantom base pair
  nframes = length(sequence);
  nbp = nframes-1;

% Call cgDNA to compute intrinsic-shape vector and stiffness matrix
  params = load('cgDNA/cgDNAparamset2.mat');
  p=path;
  path(p,'cgDNA');
  [whats, stiff] = constructSeqParms(sequence, params);
  path(p)

% Parse stiffness matrix into nonzero sub-blocks
  mmats=zeros(nbp,6,6);nmats=zeros(nframes,6,6); % (z,z) and (y,y)
  omats=zeros(nbp,6,6);pmats=zeros(nbp,6,6); % (z,y)
  qmats=zeros(nbp,6,6); % (y,neighboring y)
  for i=1:nbp
    for j=1:6
        for k=1:6
            mmats(i,j,k)=stiff(6+12*(i-1)+j, 6+12*(i-1)+k);
            nmats(i,j,k)=stiff(  12*(i-1)+j,   12*(i-1)+k);
            omats(i,j,k)=stiff(6+12*(i-1)+j,12+12*(i-1)+k);
            pmats(i,j,k)=stiff(  12*(i-1)+j, 6+12*(i-1)+k);
            qmats(i,j,k)=stiff(  12*(i-1)+j,12+12*(i-1)+k);
        end
    end
  end
  for j=1:6
    for k=1:6
        nmats(nframes,j,k)=stiff(12*nbp+j,12*nbp+k);
    end
  end
% Compute eigs of stiffness matrix
% disp('computing eigenvalues of stiffness matrix')
  stiff_mat = compute_stiffness_mat(nframes);
  small_eigs=eigs(stiff_mat,3*nframes+3*nbp,'SM');
  large_eigs=eigs(stiff_mat,3*nframes+3*nbp,'LM');
  large_eigs=large_eigs(3*nframes+3*nbp:-1:1); 
  all_log_eigs_stiff_mat=[log10(small_eigs);log10(large_eigs)];

  theta3hat = 2*pi/10.5; lk0 = round(nbp*theta3hat/(2*pi));
  zeta3hat = 3.4;
  dt = 1/nbp; R = (zeta3hat*nbp)/(2*pi);

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
  [os,qs,intras,zvec] = build_twisted_circle(nbp,tht,w);
  q4_at_1 = qs(4,nframes);

  fprintf("\nSearching for energy minimizer, %d bp, %d init Lk\n",nbp,lk);
  [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
  lk_comp = compute_lk_new(nbp,zeq);
  max_normq_err = compute_max_normq_err(nbp,zeq);
  fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

  [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
  fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));
     
  [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,all_log_eigs_stiff_mat);
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
  [os,qs,intras,zvec] = build_twisted_circle(nbp,tht,w);
  q4_at_1 = qs(4,nframes);

  fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
  [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
  lk_comp = compute_lk_new(nbp,zeq);
  max_normq_err = compute_max_normq_err(nbp,zeq);
  fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

  [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
  fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));
     
  [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,all_log_eigs_stiff_mat);
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
    [os,qs,intras,zvec] = build_twisted_circle(nbp,tht,w);
    q4_at_1 = qs(4,nframes);

    fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
    [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
    lk_comp = compute_lk_new(nbp,zeq);
    max_normq_err = compute_max_normq_err(nbp,zeq);
    fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

    [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
    fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));
     
    [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,all_log_eigs_stiff_mat);
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
      [os,qs,intras,zvec] = build_twisted_circle(nbp,tht,w);
      q4_at_1 = qs(4,nframes);

      fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
      [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
      lk_comp = compute_lk_new(nbp,zeq);
      max_normq_err = compute_max_normq_err(nbp,zeq);
      fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

      [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
      fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));
     
      [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,all_log_eigs_stiff_mat);
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
    [os,qs,intras,zvec] = build_twisted_circle(nbp,tht,w);
    q4_at_1 = qs(4,nframes);

    fprintf("\nSearching for energy minimizer, %d init Lk\n",lk);
    [zeq,temp1,temp2,all_time] = find_minimum_homegrown_returntime(nbp,zvec,[],1);
    lk_comp = compute_lk_new(nbp,zeq);
    max_normq_err = compute_max_normq_err(nbp,zeq);
    fprintf("  Link is %4.1f, |q| differs from 1 by at most %12.8f\n",lk_comp,max_normq_err)

    [eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
    fprintf("  Energy is %15.10f, normgrad is %15.10f\n",eneq,norm(gradeq));
     
    disp('computing det that contribute to J')
    [eneq,log10_overall_det,jfac_cont] = post_process(nbp,zeq,all_log_eigs_stiff_mat);
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

fclose(fid_out)