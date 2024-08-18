function Lk = compute_lk(nbp,zvec)
%
% Given a vector of unknowns "zvec" for "nbp" base pairs, compute Lk
%
% Coded by Matthew Scharf June 2017. 
% 
% Used the mathematical theory described in 
%"Computation of Writhe in Modeling of Supercoiled DNA" 
%by Konstantin Klenin Jorg Langowski, particularly equation 13 (modified 
%to calculate Linked Number instead of writhe)

%Initializing the matrices we will need

q = zeros(4,nbp);   
q(4,1) = 1;
o = zeros(3,nbp);
d1 = zeros(3,nbp);

%S1 and S2 represent the 2 DNA strands 
S1 = zeros(3,nbp+1);
S2 = zeros(3,nbp+1);

%Finding o and q by extracting zeq
for i=2:nbp
    o(:,i)=zvec(7*(i-2)+1:7*(i-2)+3);
    q(:,i)=zvec(7*(i-2)+4:7*(i-2)+7);
end

%Finding the distances of S1 and S2 from o and storing in d1
for i=1:nbp
    d1(:,i) = [q(4,i)^2+q(1,i)^2-q(2,i)^2-q(3,i)^2;
    2*q(1,i)*q(2,i)+2*q(3,i)*q(4,i);
    2*q(1,i)*q(3,i)-2*q(2,i)*q(4,i);
    ];
end

%Finding S1 and S2 using o and d1
for i=2:nbp
    S1(:,i)=o(:,i)+d1(:,i);
    S2(:,i)=o(:,i)-d1(:,i);
end

%Inserting the location of the first and last basepair into S1 and S2
S1(:,1)=d1(:,1);
S1(:,end)=S1(:,1);
S2(:,1)=-d1(:,1);
S2(:,end)=S2(:,1);

%Initializing the Linked Number
Lk=0;

%This begins our calculation of Linked Number. We modify equation 13 by
%letting i run from 1 to N instead of 2 to N and letting j run for all
%pairs on S2 instead of restricting j<i. By changing equation 13 to be
%complete in this way, we get a formula for Linked Number instead of writhe
for i=1:nbp
    
%p1 and p2 are an adjacent pair of points along S1 
    p1=S1(:,i);
    p2=S1(:,i+1);
        for j=1:nbp
            
%p3 and p4 are every adjacent pair of points along S2 which we are finding
%for each adjacent pair of points along S1
            p3=S2(:,j);
            p4=S2(:,j+1);
            
%rij represents the vector from pi to pj
            r13=p3-p1;
            r14=p4-p1;
            r23=p3-p2;
            r24=p4-p2;
            r12=p2-p1;
            r34=p4-p3;

%Finding the unit vectors orthogonal to pairings of vectors between the pairings of
%points on S1 and S2 as described by Klenin in equation 15
            n1=cross(r13,r14)/norm(cross(r13,r14));
            n2=cross(r14,r24)/norm(cross(r14,r24));
            n3=cross(r24,r23)/norm(cross(r24,r23));
            n4=cross(r23,r13)/norm(cross(r23,r13));

%Finding the omega using equation 16a
            omega=asin(dot(n1,n2))+asin(dot(n2,n3))+asin(dot(n3,n4))+asin(dot(n4,n1));
            
%Finding the sign of the omega using equation 16b
            omega=omega*sign(dot(cross(r34,r12),r13));

%Double summation of this omega using the nested "for" loops
%for each adjacent pairing along S2 and then for each 
%adjacent pairing along S1
            Lk=Lk+omega;
        end
end

%Divide by 4pi instead of 2pi as described in equation 13 in order to find
%Linked Number instead of writhe
Lk=abs(Lk/(4*pi));
return
end

%{
hold on
plot3(S1(1,:),S1(2,:),S1(3,:))
plot3(S2(1,:),S2(2,:),S2(3,:))
%}