function mat = convert_to_hess_wrt_k(zvec,hess_wrt_q)
% Convert (sparse) hessian wrt intra/o/q to (sparse) hessian wrt intra/o/k
%
% Ordering of zeq is intra1, 
%                then intra2,oq2,intra3,oq3,...oqn-2,intran-1,oqn-1
%                then intran
s = size(zvec);
zlen = s(1);
nbp = (zlen+14)/13;
newzlen=6*nbp+6*(nbp-2);
sk=10;

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp);
o = zeros(3,nbp);
for i=2:nbp-1
    o(:,i)=zvec(12+13*(i-2)+1:12+13*(i-2)+3);
    q(:,i)=zvec(12+13*(i-2)+4:12+13*(i-2)+7);
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
nintraintra = 36*nbp+36*(nbp-1)+36*(nbp-1);  % (a,a)  (a,a-1)  (a,a+1)
nintrainter = 36*(nbp-2)+36*(nbp-2)+36*(nbp-2);  % (a,a)  (a,a-1)  (a,a+1)
ninterintra = nintrainter;
ninterinter = 36*(nbp-2)+36*(nbp-3)+36*(nbp-3); % (a,a) (a,a-1) (a,a+1)
nentries = nintraintra+nintrainter+ninterintra+ninterinter;
I = zeros(nentries,1); J = zeros(nentries,1); X = zeros(nentries,1);

nset = 0;   % tracks which nonzero entry we are inputting

% (intra,intra) entries of hessian are unchanged
for a=1:nbp  % intra (a,a) slots
    if a==1
        rowrange=1:6; colrange=rowrange;
        newrowrange=1:6; newcolrange=newrowrange;
    else
        rowrange=7+13*(a-2):12+13*(a-2);colrange=rowrange;
        newrowrange=7+12*(a-2):12+12*(a-2);newcolrange=newrowrange;
    end
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=hess_wrt_q(rowrange(j),colrange(k));
        end
    end
end
for a=1:nbp-1  % intra (a,a+1) slots (and (a,a-1) done by symmetry)
    if a==1
        rowrange=1:6;colrange=rowrange+6;
        newrowrange=1:6;newcolrange=newrowrange+6;
    else
        rowrange=7+13*(a-2):12+13*(a-2); colrange=rowrange+13;
        newrowrange=7+12*(a-2):12+12*(a-2);newcolrange=newrowrange+12;
    end
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=hess_wrt_q(rowrange(j),colrange(k));
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

scal=1.0;

mat1=zeros(6,7); 
% (intra,inter) and (inter,intra) entries of hessian
for a=2:nbp-1 % (intra,inter) (a,a)   and (inter,intra) (a,a) by symmetry
    rowrange=7+13*(a-2):12+13*(a-2);colrange=rowrange(1)+6:rowrange(1)+12;
    newrowrange=7+12*(a-2):12+12*(a-2);newcolrange=newrowrange(1)+6:newrowrange(1)+11;
% Read in old Hessian
    for j=1:6
        for k=1:7
            mat1(j,k)=hess_wrt_q(rowrange(j),colrange(k));
        end
    end
% Build 7x6 matrix to right-multiply by
    mat2 = [eye(6); zeros(1,6)];
    for j=4:7
        for k=4:6
            mat2(j,k)=bq(k-3,j-3,a)/sk;
        end
    end
    result=mat1*mat2;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=result(j,k);
            if k>3
                X(nset)=X(nset)/scal;
            end
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

for a=1:nbp-2 % (intra,inter) (a,a+1)  and (inter,intra) (a+1,a) by symm
    if a==1
        rowrange=1:6;colrange=rowrange(1)+12:rowrange(1)+18;
        newrowrange=1:6;newcolrange=newrowrange(1)+12:newrowrange(1)+17;
    else
        rowrange=7+13*(a-2):12+13*(a-2);colrange=rowrange(1)+19:rowrange(1)+25;
        newrowrange=7+12*(a-2):12+12*(a-2);newcolrange=newrowrange(1)+18:newrowrange(1)+23;
    end
    % Read in old Hessian
    for j=1:6
        for k=1:7
            mat1(j,k)=hess_wrt_q(rowrange(j),colrange(k));
        end
    end
    % Build 7x6 matrix to right-multiply by
    mat2 = [eye(6); zeros(1,6)];
    for j=4:7
        for k=4:6
            mat2(j,k)=bq(k-3,j-3,a+1)/sk;
        end
    end
    result=mat1*mat2;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=result(j,k);
            if k>3
                X(nset)=X(nset)/scal;
            end
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

for a=3:nbp % (intra,inter) (a,a-1) and (inter,intra) (a-1,a) by symm
    rowrange=7+13*(a-2):12+13*(a-2);colrange=rowrange(1)-7:rowrange(1)-1;
    newrowrange=7+12*(a-2):12+12*(a-2);newcolrange=newrowrange(1)-6:newrowrange(1)-1;
    % Read in old Hessian
    for j=1:6
        for k=1:7
            mat1(j,k)=hess_wrt_q(rowrange(j),colrange(k));
        end
    end
    % Build 7x6 matrix to righbt-multiply by
    mat2 = [eye(6); zeros(1,6)];
    for j=4:7
        for k=4:6
            mat2(j,k)=bq(k-3,j-3,a-1)/sk;
        end
    end
    result=mat1*mat2;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=newrowrange(j);J(nset)=newcolrange(k);X(nset)=result(j,k);
            if k>3
                X(nset)=X(nset)/scal;
            end
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

% (inter,inter) entries of hessian
mat1=zeros(7,7); 
for a=2:nbp-2 % (a,a+1)   and (a+1,a) by symmetry
    rowrange=13+13*(a-2):19+13*(a-2);colrange=rowrange(1)+13:rowrange(1)+19;
    newrowrange=13+12*(a-2):18+12*(a-2);newcolrange=newrowrange(1)+12:newrowrange(1)+17;
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
    rowrange=13+13*(a-2):19+13*(a-2);colrange=rowrange;
    newrowrange=13+12*(a-2):18+12*(a-2);newcolrange=newrowrange;
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
mat = sparse(I,J,X,newzlen,newzlen); 