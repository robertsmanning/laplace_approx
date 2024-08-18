function energy = energy_of_twisted_circle(nbp,tht,lk)
% Given # of base pairs "nbp", register angle "tht", linking number "lk",
%   compute energy of twisted circle
  global whats
  global q4_at_1
  nframes = nbp+1; dt = 1/nbp; 
  R = (3.4*nbp)/(2*pi); w = 2*pi*lk;
  os = zeros(3,nframes); qs = zeros(4,nframes); intras = zeros(6,nframes);
  for i = 1:nframes
    os(:,i) = [R*(cos(2*pi*(i-1)*dt)-1)*sin(tht); 
            R*(cos(2*pi*(i-1)*dt)-1)*cos(tht); R*sin(2*pi*(i-1)*dt)];
    s1 = sin(pi*(i-1)*dt); c1 = cos(pi*(i-1)*dt); 
    sw = sin(w*(i-1)*dt/2); cw = cos(w*(i-1)*dt/2);
    qs(:,i) = [cos(tht)*s1*cw-sin(tht)*s1*sw;
        -sin(tht)*s1*cw-cos(tht)*s1*sw;c1*sw;c1*cw];
    intras(:,i) = whats(12*(i-1)+1:12*(i-1)+6);
    end
    zvec=zeros(6*nframes+7*(nbp-1),1);
    zvec(1:6)=intras(:,1);
    for i=2:nbp
      zvec(6+13*(i-2)+1:6+13*(i-2)+13)=[intras(:,i); os(:,i); qs(:,i)];
    end
    zvec(6+13*(nbp-1)+1:6+13*(nbp-1)+6)=intras(:,nframes);
    q4_at_1 = qs(4,nframes);
    energy = discrete_dna_penalty_en(zvec);
    return
end