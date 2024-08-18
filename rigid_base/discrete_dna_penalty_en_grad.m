function [energy,grad] = discrete_dna_penalty_en_grad(zvec)
% Given vector of unknowns "zvec", compute energy and gradient
global whats
global mmats
global nmats
global omats
global pmats
global qmats
global q4_at_1
% We have 6(nbp) intra unknowns and 7(nbp-2) unknowns (o,q),
%    with the first (o,q)=(0,e4) and the last (0,\pm e4)
%
% Ordering of zvec is intra1, 
%                then intra2,oq2,intra3,oq3,...oqn-2,intran-1,oqn-1
%                then intran
s = size(zvec);
zlen = s(1);
nbp = (zlen+14)/13;
penalty_weight=100;   % cofficient for the n-2 penalties on (|q|^2-1)^2

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp);
o = zeros(3,nbp);
for i=2:nbp-1
    o(:,i)=zvec(12+13*(i-2)+1:12+13*(i-2)+3);
    q(:,i)=zvec(12+13*(i-2)+4:12+13*(i-2)+7);
end
q(4,1)=1;
q(4,nbp)=q4_at_1;

% Read intra values from zvec
intras = zeros(6,nbp); intradiffs = zeros(6,nbp);
for i=2:nbp-1
    intradiffs(:,i)=zvec(6+13*(i-2)+1:6+13*(i-2)+6);
end
intradiffs(:,1)=zvec(1:6);
intradiffs(:,nbp)=zvec(zlen-5:zlen);
for i=1:nbp
    intras(:,i)=intradiffs(:,i)+whats(12*(i-1)+1:12*(i-1)+6);
end

% Set up key constant matrices
id = eye(4);
b = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3) = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1) = -1; b(3,3,4) = 1; b(3,4,3) = -1;

% Store values of B1q,B2q,B3q in triply-indexed array bq
bq = zeros(3,4,nbp);
for i=1:3
    bi = zeros(4,4);
    for j1 = 1:4
        for j2 = 1:4
            bi(j1,j2)=b(i,j1,j2);
        end
    end
    for k = 1:nbp
        dum = bi*q(:,k);
        for j1=1:4
            bq(i,j1,k)=dum(j1);
        end
    end
end

% Compute the inters
cay = zeros(3,nbp-1);
tr = zeros(3,nbp-1);
inters = zeros(6,nbp-1);
dfac = zeros(1,nbp-1);
for i=1:nbp-1
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bkq = zeros(4,1);
        for j1 = 1:4
            bkq(j1)=bq(k,j1,i);
        end
        cay(k,i) = 10/dfac(i)*q(:,i+1)'*bkq;
        inters(k,i)=cay(k,i);
    end
    do = o(:,i+1)-o(:,i);
    dirs = compute_ds(q(:,i+1)+q(:,i));
    for k=1:3
        tr(k,i)=do'*dirs(:,k);
        inters(k+3,i)=tr(k,i);
    end
end
% Store values of inter-interhat
interdiffs = zeros(6,nbp);
for i=1:nbp-1
    interdiffs(:,i)=inters(:,i)-whats(6+12*(i-1)+1:6+12*(i-1)+6);
end

% Compute the energy -- M contribs, then N, then O, then P, then Q
%   (no factor of 1/2 for O,P,Q since contribs from above and below diag)
energy = 0.0;
bend_energy = 0.0; twist_energy = 0.0;
mat=zeros(6,6);vec=zeros(6,1);vec2=zeros(6,1);
for i=1:nbp-1
    for j1=1:6
        vec(j1)=interdiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=mmats(i,j1,j2);
        end
    end
    energy=energy+vec'*(mat*vec)/2;
    bend_energy=bend_energy+vec(1)*mat(1,1)*vec(1)/2+vec(2)*mat(2,2)*vec(2)/2;
    twist_energy=twist_energy+vec(3)*mat(3,3)*vec(3)/2;
end
for i=1:nbp
    for j1=1:6
        vec(j1)=intradiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=nmats(i,j1,j2);
        end
    end
    energy=energy+vec'*(mat*vec)/2;
end
for i=1:nbp-1
    for j1=1:6
        vec(j1)=interdiffs(j1,i);
        vec2(j1)=intradiffs(j1,i+1);
        for j2=1:6
            mat(j1,j2)=omats(i,j1,j2);
        end
    end
    energy=energy+vec'*(mat*vec2);
end
for i=1:nbp-1
    for j1=1:6
        vec(j1)=intradiffs(j1,i);
        vec2(j1)=interdiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=pmats(i,j1,j2);
        end
    end
    energy=energy+vec'*(mat*vec2);
end
for i=1:nbp-1
    for j1=1:6
        vec(j1)=intradiffs(j1,i);
        vec2(j1)=intradiffs(j1,i+1);
        for j2=1:6
            mat(j1,j2)=qmats(i,j1,j2);
        end
    end
    energy=energy+vec'*(mat*vec2);
end


% Add penalty factors for each |q|^2
for i=1:nbp-2
    energy = energy+penalty_weight*(q(:,i+1)'*q(:,i+1)-1)^2;
end

% Fill in gradient vector
grad = zeros(zlen,1);

% Intra slots
for i=1:nbp
    contrib=zeros(6,1);
    for j1=1:6
        vec(j1)=intradiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=nmats(i,j1,j2);
        end
    end
    contrib=contrib+mat*vec;
    if i<nbp
      for j1=1:6
        vec(j1)=intradiffs(j1,i+1);
        for j2=1:6
            mat(j1,j2)=qmats(i,j1,j2);
        end
      end
      contrib=contrib+mat*vec;
    end
    if i>1
      for j1=1:6
        vec(j1)=intradiffs(j1,i-1);
        for j2=1:6
            mat(j1,j2)=qmats(i-1,j1,j2);
        end
      end
      contrib=contrib+transpose(mat)*vec;
    end
    if i<nbp
      for j1=1:6
        vec(j1)=interdiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=pmats(i,j1,j2);
        end
      end
      contrib=contrib+mat*vec;
    end
    if i>1
      for j1=1:6
        vec(j1)=interdiffs(j1,i-1);
        for j2=1:6
            mat(j1,j2)=omats(i-1,j1,j2);
        end
      end
      contrib=contrib+transpose(mat)*vec;
    end
    if i<2
        grad(1:6)=contrib;
    else
        grad(6+13*(i-2)+1:6+13*(i-2)+6)=contrib;
    end
end

% Vectors that contribute to inter slots of grad
dedz = zeros(6,nbp-1);

for i=1:nbp-1
    contrib=zeros(6,1);
    for j1=1:6
        vec(j1)=intradiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=pmats(i,j1,j2);
        end
    end
    contrib=contrib+transpose(mat)*vec;
    for j1=1:6
        vec(j1)=interdiffs(j1,i);
        for j2=1:6
            mat(j1,j2)=mmats(i,j1,j2);
        end
    end
    contrib=contrib+mat*vec;
    for j1=1:6
        vec(j1)=intradiffs(j1,i+1);
        for j2=1:6
            mat(j1,j2)=omats(i,j1,j2);
        end
    end
    contrib=contrib+mat*vec;
    dedz(:,i)=contrib;
end

% Matrices d(z^a)/d(o,q)^a that contribute to inter slots of grad
%   Note that a=1 is left as all zeroes; should not be used
dzdoq = zeros(6,7,nbp-1);

for a=2:nbp-1
    for j1=1:3
        for j2=1:4
          dzdoq(j1,j2+3,a)=(-10*bq(j1,j2,a+1)-cay(j1,a)*q(j2,a+1))/dfac(a);
        end
    end
    dirs=compute_ds(q(:,a)+q(:,a+1));
    ddirs=compute_dds(q(:,a)+q(:,a+1));
    odiff=o(:,a+1)-o(:,a);
    for j1=1:3
        for j2=1:3
            dzdoq(j1+3,j2,a)=-dirs(j2,j1);
        end
    end
    for j1=1:3
        for j2=1:4
            dzdoq(j1+3,j2+3,a)=ddirs(j1,1,j2)*odiff(1)+ddirs(j1,2,j2)*odiff(2)+ddirs(j1,3,j2)*odiff(3);
        end
    end
end
% Matrices d(z^(a-1))/d(o,q)^a that contribute to inter slots of grad
%   Note that a=1 is left as all zeroes; should not be used
dzprevdoq = zeros(6,7,nbp-1);

for a=2:nbp-1
    for j1=1:3
        for j2=1:4
          dzprevdoq(j1,j2+3,a)=(10*bq(j1,j2,a-1)-cay(j1,a-1)*q(j2,a-1))/dfac(a-1);
        end
    end
    dirs=compute_ds(q(:,a-1)+q(:,a));
    ddirs=compute_dds(q(:,a-1)+q(:,a));
    odiff=o(:,a)-o(:,a-1);
    for j1=1:3
        for j2=1:3
            dzprevdoq(j1+3,j2,a)=dirs(j2,j1);
        end
    end
    for j1=1:3
        for j2=1:4
            dzprevdoq(j1+3,j2+3,a)=ddirs(j1,1,j2)*odiff(1)+ddirs(j1,2,j2)*odiff(2)+ddirs(j1,3,j2)*odiff(3);
        end
    end
end

% Inter slots
mat67=zeros(6,7);
for a=2:nbp-1   % The index of (o,q)
    for j1=1:6
        vec(j1)=dedz(j1,a-1);
        for j2=1:7
            mat67(j1,j2)=dzprevdoq(j1,j2,a);
        end
    end
    contrib=transpose(mat67)*vec;
    for j1=1:6
        vec(j1)=dedz(j1,a);
        for j2=1:7
            mat67(j1,j2)=dzdoq(j1,j2,a);
        end
    end
    contrib=contrib+transpose(mat67)*vec;
% Add term coming from penalty weight    
    qthis=q(:,a);fac=4*penalty_weight*(qthis'*qthis-1);
    contrib=contrib+fac*[0;0;0;qthis];
    grad(12+13*(a-2)+1:12+13*(a-2)+7)=contrib;
end
