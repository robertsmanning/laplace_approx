function max_normq_err = compute_max_normq_err(nbp,zvec)
% Given a number of base pairs and a vector of unknowns, find the
%   largest error in the norm of the quaternions (each should be norm 1)
max_normq_err = 0.0;
for i=2:nbp
    q=zvec(7*(i-2)+4:7*(i-2)+7);
    normq_err = abs(norm(q)-1);
    if normq_err > max_normq_err
      max_normq_err = normq_err;
    end
end
return

