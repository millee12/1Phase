na=20;
nb=35;
state_a.type='solid';

state_a.inside.perimeter=20;
state_a.outside.perimeter=20;
state_a.mesh=rand*linspace(rand,100*rand,na+1);
state_a.N=na;
state_a.temperature=100*rand(na,1)+460;
state_a.heataddition=zeros(na,1);
state_a.R=10*ones(na,1);

state_b.type='1phasefluid';
state_b.mesh=rand*linspace(rand,100*rand,nb+1);
state_b.N=nb;
state_b.temperature=100*rand(nb,1)+460;
state_b.heataddition=zeros(nb,1);
state_b.R=1*ones(nb,1);
[lint]=calcinter3(state_a,state_b,'inside');
C=jet;
C(1,:)=1;
figure(2);
set(gca,'Fontsize',20,'LineWidth',2)
imagesc(abs(lint))
colormap(C)
h=colorbar;
set(h,'fontsize',14);
xlabel('material b')
ylabel('material a')