function result = compute_sum_log_gamma(zvec)
% Given sum of unknowns, compute sum of log of gamma values
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

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp);
for i=2:nbp-1
    q(:,i)=zvec(12+13*(i-2)+4:12+13*(i-2)+7);
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

result=0;
cay = zeros(3,nbp-1);
dfac = zeros(1,nbp-1);
for i=1:nbp-1
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bkq = zeros(4,1);
        for j1 = 1:4
            bkq(j1)=bq(k,j1,i);
        end
        cay(k,i) = 10/dfac(i)*q(:,i+1)'*bkq;
    end
    gam=1+(cay(1,i)/10)^2+(cay(2,i)/10)^2+(cay(3,i)/10)^2;
    result=result+log10(gam);
end