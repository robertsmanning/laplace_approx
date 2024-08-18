function mat = compute_stiffness_mat(nbp)
global whats
global mmats
global nmats
global omats
global pmats
global qmats
global q4_at_1
    
%   hess=0;   % Not computing this yet


% Setup for building a sparse stiffness atrix (nentries = number of nonzeros 
%    we'll put in the matrix)
nintraintra = 36*nbp+36*(nbp-1)+36*(nbp-1);  % (a,a)  (a,a-1)  (a,a+1)
nintrainter = 36*(nbp-1)+36*(nbp-1);  % (a,a)  (a,a-1)  
ninterintra = nintrainter;
ninterinter = 36*(nbp-1); % (a,a)
nentries = nintraintra+nintrainter+ninterintra+ninterinter;
I = zeros(nentries,1); J = zeros(nentries,1); X = zeros(nentries,1);

nset = 0;   % tracks which nonzero entry we are inputting

for a=1:nbp  % intra (a,a) slots
    rowrange=1+12*(a-1):6+12*(a-1);colrange=rowrange;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k); X(nset)=nmats(a,j,k);
        end
    end
end
for a=1:nbp-1  % intra (a,a+1) slots (and (a,a-1) done by symmetry)
    rowrange=1+12*(a-1):6+12*(a-1);colrange=rowrange+12;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k);X(nset)=qmats(a,j,k);
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

for a=1:nbp-1 % (intra,inter) (a,a)   and (inter,intra) (a,a) by symmetry
    rowrange=1+12*(a-1):6+12*(a-1);colrange=rowrange+6;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k);X(nset)=pmats(a,j,k);
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end
 
for a=2:nbp % (intra,inter) (a,a-1) and (inter,intra) (a-1,a) by symm
    rowrange=1+12*(a-1):6+12*(a-1);colrange=rowrange-6;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k);X(nset)=omats(a-1,k,j);
            nset=nset+1;
            I(nset)=J(nset-1);J(nset)=I(nset-1);X(nset)=X(nset-1);
        end
    end
end

for a=1:nbp-1 % (inter,inter) (a,a)
    rowrange=7+12*(a-1):12+12*(a-1);colrange=rowrange;
    for j=1:6
        for k=1:6
            nset=nset+1;
            I(nset)=rowrange(j);J(nset)=colrange(k); X(nset)=mmats(a,j,k);
        end
    end
end
mat = sparse(I,J,X,6*nbp+6*(nbp-1),6*nbp+6*(nbp-1)); 