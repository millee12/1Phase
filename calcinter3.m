function [Qdot]=HX2(state_a,state_b)
%Define Parameters
nela=state_a.N;
nelb=state_b.N;
Mrs=zeros(nela,nelb);
Mls=zeros(nela,nelb);
DT=zeros(nela,nelb);
U=zeros(nela,nelb);
%perimeter data must come from solid
if strcmp(state_a.type,'solid')
    per=eval(['state_a.',location,'.perimeter']);
elseif strcmp(state_b.type,'solid')
    per=eval(['state_b.',location,'.perimeter']);
end
%Calculate total interface length
Lint=min(state_a.mesh(nela+1),state_b.mesh(nelb+1))...
    -max(state_a.mesh(1),state_b.mesh(1));
%Calculate Matrix of possible interfaces
for a=1:nela
    for b=1:nelb
        Mrs(a,b)=min(state_a.mesh(a+1),state_b.mesh(b+1));
        Mls(a,b)=max(state_a.mesh(a),state_b.mesh(b));
        DT(a,b)=state_a.temperature(a)-state_b.temperature(b);
        U(a,b)=1/(state_a.R(a)+state_b.R(b));
    end
end
Mint=Mrs-Mls;
%Remove Negative interface Lengths
Mint(Mint<0)=0;
err=abs(sum(sum(Mint))-Lint)/Lint;
%Check Interface
if err > .1
    disp('warning: large residual in interface calculator:')
    disp(['residual= ',num2str(err)])
end
Qdot=per*Mint.*U.*DT;