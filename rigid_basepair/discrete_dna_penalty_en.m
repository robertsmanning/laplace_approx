function energy = discrete_dna_penalty_en(zvec)
% Given vector "zvec" of all unknowns, compute its energy
global stiffmat
global whats
global q4_at_1
% We have 7(nbp-2) unknowns (o,q),
% Ordering of zvec is oq2,oq3,...oqn-1
s = size(zvec);
zlen = s(1);
nbp = (zlen+14)/7;
penalty_weight=100;   % cofficient for the n-1 penalties on (|q|^2-1)^2

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp);
o = zeros(3,nbp);
for i=2:nbp-1
    o(:,i)=zvec(7*(i-2)+1:7*(i-2)+3);
    q(:,i)=zvec(7*(i-2)+4:7*(i-2)+7);
end
q(4,1)=1;
q(4,nbp)=q4_at_1;

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
cay = zeros(3,nbp-1); tr = zeros(3,nbp-1);
inters = zeros(6,nbp-1); dfac = zeros(1,nbp-1);
for i=1:nbp-1
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bkq = zeros(4,1);
        for j1 = 1:4
            bkq(j1)=bq(k,j1,i);
        end
        cay(k,i) = 2/dfac(i)*q(:,i+1)'*bkq;
        inters(k,i)=cay(k,i);
    end
    do = o(:,i+1)-o(:,i);
    dirs = compute_ds(q(:,i+1)+q(:,i));
    for k=1:3
        tr(k,i)=do'*dirs(:,k);
        inters(k+3,i)=tr(k,i);
    end
end

% Compute the energy
energy = 0.0;
for i=1:nbp-1
    vec = inters(:,i)-whats(:,i);
    energy = energy + 0.5*vec'*stiffmat*vec;
    energy = energy+penalty_weight*(q(:,i+1)'*q(:,i+1)-1)^2;
end
