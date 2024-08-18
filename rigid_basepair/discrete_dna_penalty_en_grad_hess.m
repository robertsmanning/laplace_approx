function [energy,grad,hess] = discrete_dna_penalty_en_grad_hess(zvec)
% Given vector "zvec" of all unknowns, compute its energy, gradient, and
% Hessian

global stiffmat
global whats
global q4_at_1
s = size(zvec);
zlen = s(1);
nbp = round(zlen/7)+1;
nframes = nbp+1;
% We have 7(nbp-1) unknowns (o,q),
%    with the first (o,q)=(0,e4) and the last (0,\pm e4)
%
% Ordering of zvec is oq2,oq3,...oqn

penalty_weight=100;   % cofficient for the n-1 penalties on (|q|^2-1)^2

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nframes);
o = zeros(3,nframes);
for i=2:nbp
    o(:,i)=zvec(7*(i-2)+1:7*(i-2)+3);
    q(:,i)=zvec(7*(i-2)+4:7*(i-2)+7);
end
q(4,1)=1;
q(4,nframes)=q4_at_1;

% Set up key constant matrices
id = eye(4);
b = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3) = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1) = -1; b(3,3,4) = 1; b(3,4,3) = -1;

% Store values of B1q,B2q,B3q in triply-indexed array bq
bq = zeros(3,4,nframes);
for i=1:3
    bi = zeros(4,4);
    for j1 = 1:4
        for j2 = 1:4
            bi(j1,j2)=b(i,j1,j2);
        end
    end
    for k = 1:nframes
        dum = bi*q(:,k);
        for j1=1:4
            bq(i,j1,k)=dum(j1);
        end
    end
end

% Compute the inters
cay = zeros(3,nbp); tr = zeros(3,nbp);
inters = zeros(6,nbp); interdiffs = zeros(6,nbp);
dfac = zeros(1,nbp);
for i=1:nbp
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
    interdiffs(:,i)=inters(:,i)-whats(:,i);
end

% Compute the energy 
energy = 0.0;
for i=1:nbp
    vec = inters(:,i)-whats(:,i);
    energy = energy + 0.5*vec'*stiffmat*vec;
    energy = energy+penalty_weight*(q(:,i+1)'*q(:,i+1)-1)^2;
end

% Fill in gradient vector
grad = zeros(zlen,1);

% Vectors that contribute to inter slots of grad (and to Hessian)
dedz = zeros(6,nbp);
for i=1:nbp
    dedz(:,i)=stiffmat*interdiffs(:,i);
end

% Matrices d(z^a)/d(o,q)^a that contribute to inter slots of grad (and to Hessian)
%   Note that a=1 is left as all zeroes; should not be used
dzdoq = zeros(6,7,nbp);

for a=2:nbp
    for j1=1:3
        for j2=1:4
          dzdoq(j1,j2+3,a)=(-2*bq(j1,j2,a+1)-cay(j1,a)*q(j2,a+1))/dfac(a);
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
% Matrices d(z^(a-1))/d(o,q)^a that contribute to inter slots of grad (and to Hessian)
%   Note that a=1 is left as all zeroes; should not be used
dzprevdoq = zeros(6,7,nbp);

for a=2:nbp
    for j1=1:3
        for j2=1:4
          dzprevdoq(j1,j2+3,a)=(2*bq(j1,j2,a-1)-cay(j1,a-1)*q(j2,a-1))/dfac(a-1);
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
for a=2:nbp   % The index of (o,q)
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
    grad(7*(a-2)+1:7*(a-2)+7)=contrib;
end
    
%   hess=0;   
% Setup for building a sparse Hessian (nentries = number of nonzeros 
%    we'll put in the Hessian)
ninterinter = 49*(nbp-1)+49*(nbp-2)+49*(nbp-2); % (a,a) (a,a-1) (a,a+1)
nentries = ninterinter;
I = zeros(nentries,1); J = zeros(nentries,1); X = zeros(nentries,1);

nset = 0;   % tracks which nonzero entry we are inputting

% Some new matrices to prepare for (inter,inter) derivs
d2thdq2=zeros(nbp-1,3,4,4);  % a=1 left as zero (will not be used)
d2thdqnext2=zeros(nbp-2,3,4,4);
d2thdqnextdq=zeros(nbp-2,3,4,4); % a=1 left as zero (will not be used)
%  NOTE: in notation of paper, should probably be called d2thdqdqnext
%   but not changing since already coded with this name
dodotd2s=zeros(nbp-1,3,4,4);
for a=1:nbp
    qthis=q(:,a);qnext=q(:,a+1);dotprod=transpose(qnext)*qthis;
    othis=o(:,a);onext=o(:,a+1);d2s=compute_drdotd2ds(qthis+qnext,onext-othis);
    for j=1:3
        caythis=cay(j,a);
        bqthis=zeros(4,1);bqnext=zeros(4,1);
        for k=1:4
            bqthis(k)=bq(j,k,a);bqnext(k)=bq(j,k,a+1);
        end
        %  (a,a) deriv
        mat=2*caythis*(qnext*qnext')+2*(bqnext*qnext')+2*(qnext*bqnext');
        mat=mat/(dotprod*dotprod);
        if a>1
            for k1=1:4
                for k2=1:4
                    d2thdq2(a,j,k1,k2)=mat(k1,k2);
                end
            end
        end
        % (a+1,a+1) deriv
        mat=2*caythis*(qthis*qthis')-2*(bqthis*qthis')-2*(qthis*bqthis');
        mat=mat/(dotprod*dotprod);
        if a<nbp
            for k1=1:4
                for k2=1:4
                    d2thdqnext2(a,j,k1,k2)=mat(k1,k2);
                end
            end
        end
        % (a+1,a) deriv
        mat=2*caythis*(qnext*qthis')-2*(qnext*bqthis')+2*(bqnext*qthis');
        mat=mat/(dotprod*dotprod);
        if a<nbp && a>1
            for k1=1:4
                for k2=1:4
                    d2thdqnextdq(a,j,k1,k2)=mat(k1,k2)+(-2*b(j,k1,k2)-caythis*id(k1,k2))/dotprod;
                end
            end
        end
        % dot products of do with second derivs of dj
        for k1=1:4
            for k2=1:4
                dodotd2s(a,j,k1,k2)=d2s(j,k1,k2);
            end
        end
    end
end

% (inter,inter) entries of hessian
mat1=zeros(7,7);mat2=zeros(7,6);mat3=zeros(6,6);mat4=zeros(6,7); 
for a=2:nbp-1 % (a,a+1)   and (a+1,a) by symmetry
    rowrange=1+7*(a-2):7+7*(a-2);colrange=rowrange(1)+7:rowrange(1)+13;
    vec=zeros(6,1);
    for j=1:6
        vec(j)=dedz(j,a);
    end
    ddirs=compute_dds(q(:,a)+q(:,a+1));
    for n=1:3
        mat=zeros(7,6);
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=ddirs(k2-3,n,k1-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    for n=4:7
        mat=zeros(7,6);
        for k1=1:3
            for k2=4:6
                mat(k1,k2)=-ddirs(k2-3,k1,n-3);
            end
        end
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=dodotd2s(a,k2-3,n-3,k1-3);
            end
        end
        for k1=4:7
            for k2=1:3
                mat(k1,k2)=d2thdqnextdq(a,k2,k1-3,n-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    mat3 = stiffmat;
    for k1=1:6
        for k2=1:7
            mat4(k1,k2)=dzprevdoq(k1,k2,a+1);
            mat2(k2,k1)=dzdoq(k1,k2,a);
        end
    end
    result=mat1+mat2*mat3*mat4;
    for j=1:7
        for k=1:7
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k);X(nset)=result(j,k);
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

for a=2:nbp % (a,a)
    rowrange=1+7*(a-2):7+7*(a-2);colrange=rowrange;
    vec=zeros(6,1);
    for j=1:6
        vec(j)=dedz(j,a-1);
    end
    ddirs=compute_dds(q(:,a-1)+q(:,a));
    for n=1:3
        mat=zeros(7,6);
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=ddirs(k2-3,n,k1-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    for n=4:7
        mat=zeros(7,6);
        for k1=1:3
            for k2=4:6
                mat(k1,k2)=ddirs(k2-3,k1,n-3);
            end
        end
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=dodotd2s(a-1,k2-3,n-3,k1-3);
            end
        end
        for k1=4:7
            for k2=1:3
                mat(k1,k2)=d2thdqnext2(a-1,k2,k1-3,n-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    mat3 = stiffmat;
    for k1=1:6
        for k2=1:7
            mat4(k1,k2)=dzprevdoq(k1,k2,a);
            mat2(k2,k1)=dzprevdoq(k1,k2,a);
        end
    end
    result=mat1+mat2*mat3*mat4; 
 % Second half of formula   
    vec=zeros(6,1);
    for j=1:6
        vec(j)=dedz(j,a);
    end
    ddirs=compute_dds(q(:,a)+q(:,a+1));
    for n=1:3
        mat=zeros(7,6);
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=-ddirs(k2-3,n,k1-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    for n=4:7
        mat=zeros(7,6);
        for k1=1:3
            for k2=4:6
                mat(k1,k2)=-ddirs(k2-3,k1,n-3);
            end
        end
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=dodotd2s(a,k2-3,n-3,k1-3);
            end
        end
        for k1=4:7
            for k2=1:3
                mat(k1,k2)=d2thdq2(a,k2,k1-3,n-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    mat3 = stiffmat;
    for k1=1:6
        for k2=1:7
            mat4(k1,k2)=dzdoq(k1,k2,a);
            mat2(k2,k1)=dzdoq(k1,k2,a);
        end
    end
    result=result+mat1+mat2*mat3*mat4; 
% Penalty weight part of formula
    qthis=q(:,a);vec=[0;0;0;qthis];
    mat5=eye(7);mat5(1,1)=0;mat5(2,2)=0;mat5(3,3)=0;mat5=4*penalty_weight*(qthis'*qthis-1)*mat5;
    mat5=mat5+8*penalty_weight*(vec*vec');
    result=result+mat5;
    for j=1:7
        for k=1:7
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k);X(nset)=result(j,k);
        end
    end 
end

hess = sparse(I,J,X,zlen,zlen);        
        
   

