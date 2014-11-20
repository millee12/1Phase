function [vinterp]= interp2DUG(tvew,tvsn,M,pvew,pvsn)
n=length(pvsn);
indmax=length(tvsn);
%Vector Gradient (assumed uniform)
dtvew=tvew(2)-tvew(1);
PEWind=((pvew-tvew(1))./dtvew)+1;
PEWind(PEWind >= indmax)=indmax-.5;
PEWind(PEWind <= 1)=1.5;
%East Parameters
PE=floor(PEWind);
%West Parameters
PW=PE+1;
%South-North Parameters
dtvsn=tvsn(2)-tvsn(1);
PSNind=((pvsn-tvsn(1))./dtvsn)+1;
PSNind(PSNind >= indmax)=indmax-.5;
PSNind(PSNind <= 1)=1.5;
%South Parameters
PS=floor(PSNind);
%North Parameters
PN=ceil(PSNind);
vinterp=zeros(n,1);
for i=1:n
    %Estimate partials
    dMSN=(M(PN(i),PE(i))-M(PS(i),PE(i)))*(PSNind(i)-PS(i));
    dMEW=(M(PS(i),PW(i))-M(PS(i),PE(i)))*(PEWind(i)-PE(i));
    %add partials to interpolated value
    vinterp(i)=M(PS(i),PE(i))+dMSN+dMEW;
end

