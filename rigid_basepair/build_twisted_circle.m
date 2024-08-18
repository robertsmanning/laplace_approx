function [os,qs,zvec] = build_twisted_circle(nbp,tht,w)
% Given # of base pairs "nbp", register angle "tht", linking number "lk",
%   build twisted circle (outputs are origins "os", quaternions "qs",
%   and vector of all unknowns "zvec")
  global whats
  nframes = nbp+1; dt = 1/nbp; R = (3.4*nbp)/(2*pi);
  os = zeros(3,nframes); qs = zeros(4,nframes);
  for i = 1:nframes
    os(:,i) = [R*(cos(2*pi*(i-1)*dt)-1)*sin(tht); 
            R*(cos(2*pi*(i-1)*dt)-1)*cos(tht); R*sin(2*pi*(i-1)*dt)];
    s1 = sin(pi*(i-1)*dt); c1 = cos(pi*(i-1)*dt); 
    sw = sin(w*(i-1)*dt/2); cw = cos(w*(i-1)*dt/2);
    qs(:,i) = [cos(tht)*s1*cw-sin(tht)*s1*sw;
            -sin(tht)*s1*cw-cos(tht)*s1*sw;c1*sw;c1*cw];
  end
  zvec=zeros(7*(nbp-1),1);
  for i=2:nbp
    zvec(7*(i-2)+1:7*(i-2)+7)=[os(:,i); qs(:,i)];
  end
  return
end
