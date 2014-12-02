function [state,Ind]=remesh(state)
%Define Variables
bp=state.geometry.breakpoints;
ns=length(bp)-1;
vf=state.mass./state.density;
nf=length(vf);
A=state.geometry.crosssections;
%Define Solid section lengths and volumes
ls=bp(2:ns+1)-bp(1:ns);
Ls=sum(ls);
vs=A.*ls;
%generate segment indicator matrix
Ind=zeros(nf,ns);
a=1;b=1;
L=zeros(nf,ns);
while a<=nf && b<=ns
    Vf=sum(vf(1:a));
    Vs=sum(vs(1:b));
    L(a,b)=vf(a)/A(b);
    res=(Vf-Vs)/vf(a);
    if res <= 0
        Ind(a,b)=1-sum(sum(Ind(a,:)));
        a=a+1;
    else
        Ind(a,b)=1-res;
        b=b+1;
    end
end
Lf=L.*Ind;
state.mesh=zeros(nf+1,1);
state.mesh(1)=bp(1);
for a=1:nf
    state.mesh(a+1)=state.mesh(a)+sum(Lf(a,:));
end
check=abs(sum(sum(Lf))-Ls)/Ls;
if check > .01
    disp(['warning, mesh residual high in ',state.name])
       disp(['length discrepancy=',num2str(check)])
end
