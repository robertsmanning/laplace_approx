function [shapes, stiff] = constructSeqParmsPrime(seq, params)

%-----------------------------------------------------------
% cgDNA function: [shapes,stiff] = constructSeqParms(seq, params)
%-----------------------------------------------------------
% This function constructs the ground-state coordinate
% vector and stiffness matrix in non-dimensional Curves+
% form for a given sequence, using the specified parameter
% set in params. Note that params have to be in the 'Prime' 
% format.
%
%
% Input: 
%
%   seq     sequence along reference strand
%
%   params  parameter set structure (see Note 1). 
%
%
% Output:
%
%   shapes  ground-state coordinate vector 
%           [size N x 1]
%
%   stiff   ground-state stiffness matrix
%           [size N x N]
%
%   where N = 12*nbp - 6 and nbp is the length 
%   of the sequence seq (number of basepairs). 
%
%
% Note 1:
%
%    'params' ('Prime' format!) is an object with the 'dimer' 
%    properties defined, containing stiffness blocks and 
%    weighted shape vectors for each of the 16 possible 
%    base-pair steps in the middle of the sequence, 
%    16 base-pair steps at 5' end and 16 base-pair steps at 3' end. 
%    
%    - 'dimer': a 1x16 struct array with fields: 
%      - 'S'  : the basepair step as a string of 2 chars;
%      - 'b18': the 18x18 stiffness block corresponding 
%               to basepair 'S',
%      - 'c18': the 18x1 weighted shape vector corresponding
%               to basepair 'S',
%      - 'end1b18': 18x18 stiffness block corresponding 
%               to a 5' end basepair 'S',
%      - 'end1c18': 18x1 weighted shape vector corresponding
%               to a 5' end basepair 'S',
%      - 'end2b18': 18x18 stiffness block corresponding 
%               to a 3' end basepair 'S',
%      - 'end2c18': 18x1 weighted shape vector corresponding 
%               to a 3' end basepair 'S',
%
% Note 2:
%
%    The entries in the variables shapes and stiff
%    are consistent with the following ordering of
%    the structural coordinates
%
%     y_1, z_1, ..., y_{nbp-1}, z_{nbp-1}, y_{nbp}
%
%    where for each a=1,2,3,... we have
%
%     y_a = (Buckle,Propeller,Opening,Shear,Stretch,Stagger)_a
%
%     z_a = (Tilt,Roll,Twist,Shift,Slide,Rise)_a.
%
%    For example
%
%     shapes((a-1)*12+1) = Buckle at basepair a
%     shapes((a-1)*12+2) = Propeller at basepair a
%      ...
%     shapes((a-1)*12+6) = Stagger at basepair a
%     shapes((a-1)*12+7) = Tilt at junction a, a+1
%     shapes((a-1)*12+8) = Roll at junction a, a+1
%      ...
%     shapes((a-1)*12+12) = Rise at junction a, a+1.
%
%    Correspondingly, we have
%
%     stiff(i,j) = stiffness coefficient for the pair of
%                  coordinates shapes(i) and shapes(j).
%
%
% If you find this code useful, please cite:
%--------------------------------------------------------

    dimer = params.dimer;

    % Initialize variables
    seq = upper(seq);
    nbp = numel(seq);
    nv = 12; nover = 6;
    N = nv*(nbp-1)+nover;
    nz = (nv^2 + 2*nv*nover) * (nbp-1) + nover^2;
    stiff = spalloc(N,N,nz);
    sigma = zeros(N,1);


    % Assemble stiffness matrix 'stiff' and sigma vector 'sigma'
    %fprintf('Constructing stiffness matrix...\n');tic;
    
    stiff(1:18,1:18) = stiff(1:18,1:18) + dimer(fsi(dimer,seq(1:2))).end1b18;
    sigma(1:18) = sigma(1:18) + dimer(fsi(dimer,seq(1:2))).end1c18;
        
    for i = 2:nbp-2
        
      k = (i-1)*12 + 1;
      stiff(k:k+17,k:k+17) = stiff(k:k+17,k:k+17) + dimer(fsi(dimer,seq(i:i+1))).b18;
      sigma(k:k+17) = sigma(k:k+17) + dimer(fsi(dimer,seq(i:i+1))).c18;
    
    end
      
    stiff(N-17:N,N-17:N) = stiff(N-17:N,N-17:N) + dimer(fsi(dimer,seq(nbp-1:nbp))).end2b18;
    sigma(N-17:N) = sigma(N-17:N) + dimer(fsi(dimer,seq(nbp-1:nbp))).end2c18;

    % Compute ground-state coord vector via matrix inversion 
    %fprintf('Constructing ground-state shape...\n');toc
    shapes = stiff\sigma;

end

%--------------------------------------------------------
function i = fsi(struc, s)

    n = size(struc);

    for j = 1:n(2)
        if(strcmpi(struc(j).S,s)) 
            i=j; 
            return;
        end;    
    end    

end
