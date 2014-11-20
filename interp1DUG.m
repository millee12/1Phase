function [vinterp]= interp1DUG(tvsn,M,pvsn)
n=length(pvsn);
indmax=length(tvsn);
%South-North Parameters
dtvsn=tvsn(2)-tvsn(1);
PSNind=((pvsn-tvsn(1))./dtvsn)+1;
PSNind(PSNind >= indmax)=indmax-.5;
PSNind(PSNind <= 1)=1.5;
%South Parameters
PS=floor(PSNind);
%North Parameters
PN=ceil(PSNind);
%Interpolated Value
vinterp=zeros(n,1);
for i=1:n
    %Estimate partials
    dMSN=(M(PN(i))-M(PS(i)))*(PSNind(i)-PS(i));
    %add partials to interpolated value
    vinterp(i)=M(PS(i))+dMSN;
end

