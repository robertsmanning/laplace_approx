function energy = discrete_dna_penalty_en(zvec)
% Given vector of unknowns "zvec", compute energy
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
%max_intradiffs=max(max(abs(intradiffs)))

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
for i=1:nbp-1
    energy = energy+penalty_weight*(q(:,i+1)'*q(:,i+1)-1)^2;
end