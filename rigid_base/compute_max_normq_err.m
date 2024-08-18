function max_normq_err = compute_max_normq_err(nbp,zvec)
% Given vector of unknowns "zvec" for "nbp" base pairs, find
%   max deviation of any quaternion norm from 1
max_normq_err = 0.0;
for i=2:nbp-1
    q=zvec(12+13*(i-2)+4:12+13*(i-2)+7);
    normq_err = abs(norm(q)-1);
    if normq_err > max_normq_err
      max_normq_err = normq_err;
    end
end
return

