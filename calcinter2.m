function [Mint]=calcinter2(mesh_a,mesh_b)
%Define Parameters
na=length(mesh_a);
nela=na-1;
nb=length(mesh_b);
nelb=nb-1;
%Calculate total interface length
Lint=min(mesh_a(na),mesh_b(nb))-max(mesh_a(1),mesh_b(1));
%Calculate Leftmost Points
lsa=mesh_a(1:na-1);
lsb=mesh_b(1:nb-1);
%Calculate Rightmost points
rsa=mesh_a(2:na);
rsb=mesh_b(2:nb);
%Calculate Matrix of possible interfaces
Ma=ones(nelb,1)*rsa';
Mb=rsb*ones(1,nela);
Mrs=min(Ma,Mb);
Ma=ones(nelb,1)*lsa';
Mb=lsb*ones(1,nela);
Mls=max(Ma,Mb);
Mint=Mrs-Mls;
%Remove Negative interface Lengths
Mint(Mint<0)=0;
err=abs(sum(sum(Mint))-Lint)/Lint;
%Check Result
if err > .1
    disp('warning: large residual in interface calculator:')
    disp(['residual= ',num2str(err)])
end