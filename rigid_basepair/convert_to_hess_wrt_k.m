function mat = convert_to_hess_wrt_k(zvec,hess_wrt_q)
% Convert (sparse) hessian wrt intra/o/q to (sparse) hessian wrt intra/o/k
%
% Ordering of zeq is intra1, 
%                then intra2,oq2,intra3,oq3,...oqn-2,intran-1,oqn-1
%                then intran
s = size(zvec);
zlen = s(1);
nbp = (zlen+14)/7;
newzlen=6*(nbp-2);
sk=2;

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp);
o = zeros(3,nbp);
for i=2:nbp-1
    o(:,i)=zvec(7*(i-2)+1:7*(i-2)+3);
    q(:,i)=zvec(7*(i-2)+4:7*(i-2)+7);
end
q(4,1)=1;
q(4,nbp)=1;   % Could be -1 but this never gets used anyway

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
    
%   hess=0;   % Not computing this yet


% Setup for building a sparse Hessian (nentries = number of nonzeros 
%    we'll put in the new Hessian)
ninterinter = 36*(nbp-2)+36*(nbp-3)+36*(nbp-3); % (a,a) (a,a-1) (a,a+1)
nentries = ninterinter;
I = zeros(nentries,1); J = zeros(nentries,1); X = zeros(nentries,1);

nset = 0;   % tracks which nonzero entry we are inputting

scal=1.0;

% (inter,inter) entries of hessian
mat1=zeros(7,7); 
for a=2:nbp-2 % (a,a+1)   and (a+1,a) by symmetry
    rowrange=1+7*(a-2):7+7*(a-2);colrange=rowrange(1)+7:rowrange(1)+13;
    newrowrange=1+6*(a-2):6+6*(a-2);newcolrange=newrowrange(1)+6:newrowrange(1)+11;
    % Read in old Hessian
    for j=1:7
        for k=1:7
            mat1(j,k)=hess_wrt_q(rowrange(j),colrange(k));
        end
    end
    % Build 7x6 matrix we right-multiply by
    mat2 = [eye(6); zeros(1,6)];
    for j=4:7
        for k=4:6
            mat2(j,k)=bq(k-3,j-3,a+1)/sk;
        end
    end
    % Build 7x6 matrix whose transpose we left-multiply by
    mat3 = [eye(6); zeros(1,6)];
    for j=4:7
        for k=4:6
            mat3(j,k)=bq(k-3,j-3,a)/sk;
        end
    end
    result=transpose(mat3)*mat1*mat2;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=result(j,k);
            if k>3
                X(nset)=X(nset)/scal;
            end
            if j>3
                X(nset)=X(nset)/scal;
            end
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

for a=2:nbp-1 % (a,a)
    rowrange=1+7*(a-2):7+7*(a-2);colrange=rowrange;
    newrowrange=1+6*(a-2):6+6*(a-2);newcolrange=newrowrange;
    % Read in old Hessian
    for j=1:7
        for k=1:7
            mat1(j,k)=hess_wrt_q(rowrange(j),colrange(k));
        end
    end
    % Build 7x6 matrix we right-multiply by, and whose transpose we
    % left-multiply by
    mat2 = [eye(6); zeros(1,6)];
    for j=4:7
        for k=4:6
            mat2(j,k)=bq(k-3,j-3,a)/sk;
        end
    end
    result=transpose(mat2)*mat1*mat2; 
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=result(j,k);
            if k>3
                X(nset)=X(nset)/scal;
            end
            if j>3
                X(nset)=X(nset)/scal;
            end               
        end
    end 
end
%nset
% zlen
% max(I)
% max(J)
mat = sparse(I,J,X,newzlen,newzlen); 