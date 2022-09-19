%% This code compares connection numbers and R synchronization in %distant-dependent and probablity-dependent model
dmax=0:0.1:2.6;
connection_d=zeros(1,length(dmax));
for i=1:length(dmax)
connection_d(i)=cal_connection_d(dmax(i));
end

p_rec=0:0.1:1; 
connection_p=zeros(1,length(p_rec));
for i=1:length(p_rec)
connection_p(i)=cal_connection_p(p_rec(i));
end

figure
plot(dmax/2.6,connection_d,'linewidth',2)
hold on
plot(p_rec,connection_p,'linewidth',2)
ylabel('Average number of links per shell neuron')
xlabel('Normalized distance threshold or probability')
legend('Distance-dependent coupling','Probability-dependent coupling')
set(gca,'FontSize',14)
grid on
%%
figure
dmax=[0 0.5 1 1.5 2 2.5 3];
p_rec=[0 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 1]; 
plot(dmax/3, R_rec_d(2,:),'-o','linewidth',2)
hold on
plot(p_rec, mean(Rs),'-o','linewidth',2)
xlabel('Normalized distance threshold or probability')
ylabel('R_s_y_n (shell)')
ylim([0 0.9])
%legend('Distance-dependent coupling','Probability-dependent coupling')
set(gca,'FontSize',14)
grid on
%%
connection_d=[0	2.905555556	19.46666667	41.80555556	55.71111111	59.93333333	60];
connection_p=[0	0.06	0.6	3	6	12	18	24	30	36	60];
figure
plot(connection_d,R_rec_d(2,:),'-o','linewidth',2)
hold on
plot(connection_p,mean(Rs),'-o','linewidth',2)
xlabel('Average number of links per shell neuron')
ylabel('R_s_y_n (shell)')
ylim([0 0.9])
%legend('Diffusible coupling network','Random coupling network')
set(gca,'FontSize',14)
grid on
%% dmax
n=60; 
R=0.5; x0=1; y0=1;
rng(n)
 t = 2*pi*rand(n,1);
 r = R*sqrt(rand(n,1));
 x = x0 + r.*cos(t);
 y = y0 + r.*sin(t);
 

n1=5*n;R1=0.75;R2=1/R1;x01=1; y01=1.5+R2-1;
rng(n1)
t = 2*pi*rand(n1,1);
r1 = R1*sqrt(rand(n1,1)); r2 = r1*R2/R1;
x1 = x01 + r1.*cos(t);
y1 = y01 + r2.*sin(t);
xy1_rec=[]; 
 for i=1:n1
 if ((x1(i)-1)^2+(y1(i)-1)^2)^0.5>R
      xy1_rec=[xy1_rec;x1(i),y1(i)];
 end
 end
n2=3*n;
location=[x,y];
location2=xy1_rec(1:3*n,:);
figure
plot(location(:,1),location(:,2),'o','Linewidth',2)
hold on
plot(location2(:,1),location2(:,2),'o','Linewidth',2)


dmax=0.3;
for i=1:n
for j=1:n2
    if ((location2(j,1)-location(i,1))^2+(location2(j,2)-location(i,2))^2)^0.5<dmax
        plot([location(i,1),location2(j,1)],[location(i,2),location2(j,2)],'--','color','k')
    end
end
end

for i=1:n2
for j=1:n2
    if ((location2(i,1)-location2(j,1))^2+(location2(i,2)-location2(j,2))^2)^0.5<0.2
        plot([location2(i,1),location2(j,1)],[location2(i,2),location2(j,2)],'--','color','r')
    end
end
end 

for i=1:n
for j=1:n
    if ((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5<0.2
        plot([location(i,1),location(j,1)],[location(i,2),location(j,2)],'--','color','b')
    end
end
end
axis off ; 

% degree distribution
shell_degree=zeros(1,n2);
%dmax=0.8; %connection=10
%dmax=1; %connection=20
dmax=1.5; %connection=40
for i=1:n2
for j=1:n
    if ((location(j,1)-location2(i,1))^2+(location(j,2)-location2(i,2))^2)^0.5<dmax
shell_degree(i)=shell_degree(i)+1;
    end
end
end
figure
histogram(shell_degree,'facecolor',[0 0 0]+0.5,'BinWidth',1,'normalization','probability');xlim([0,60]);ylim([0,0.5])
xlabel('Average number of links per shell neuron')
ylabel('Normalized frequency')
set(gca,'FontSize',14)
%% p
n=60; 
R=0.5; x0=1; y0=1;
rng(n)
 t = 2*pi*rand(n,1);
 r = R*sqrt(rand(n,1));
 x = x0 + r.*cos(t);
 y = y0 + r.*sin(t);
 

n1=5*n;R1=0.75;R2=1/R1;x01=1; y01=1.5+R2-1;
rng(n1)
t = 2*pi*rand(n1,1);
r1 = R1*sqrt(rand(n1,1)); r2 = r1*R2/R1;
x1 = x01 + r1.*cos(t);
y1 = y01 + r2.*sin(t);
xy1_rec=[]; 
 for i=1:n1
 if ((x1(i)-1)^2+(y1(i)-1)^2)^0.5>R
      xy1_rec=[xy1_rec;x1(i),y1(i)];
 end
 end
n2=3*n;
location=[x,y];
location2=xy1_rec(1:3*n,:);
figure
plot(location(:,1),location(:,2),'o','Linewidth',2)
hold on
plot(location2(:,1),location2(:,2),'o','Linewidth',2)
p=0.69/60;
p=0.2/60;

for i=1:n2
for j=1:n2
    if ((location2(i,1)-location2(j,1))^2+(location2(i,2)-location2(j,2))^2)^0.5<0.2
        plot([location2(i,1),location2(j,1)],[location2(i,2),location2(j,2)],'--','color','r')
    end
end
end

for i=1:n
for j=1:n
    if ((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5<0.2
        plot([location(i,1),location(j,1)],[location(i,2),location(j,2)],'--','color','b')
    end
end
end
 rng(12)
for i=1:n
for j=1:n2
   xx=rand;
if xx<p
plot([location(i,1),location2(j,1)],[location(i,2),location2(j,2)],'-','color','k','LineWidth',2)
end
end
end

axis off ; 

% degree distribution

%p=10/60; %connection=10
p=20/60; %connection=20
%p=40/60; %connection=40
shell_degree=zeros(1,n2);
for i=1:n2
for j=1:n
   xx=rand;
if xx<p
shell_degree(i)=shell_degree(i)+1;
end
end
end

figure
histogram(shell_degree,'facecolor',[0 0 0]+0.5,'BinWidth',1,'normalization','probability');xlim([0,60]);ylim([0,0.2])
xlabel('Average number of links per shell neuron')
ylabel('Normalized frequency')
set(gca,'FontSize',14)
%% dmax
function connection_d_ave=cal_connection_d(dmax)
n=60; 
R=0.5; x0=1; y0=1;
rng(n)
 t = 2*pi*rand(n,1);
 r = R*sqrt(rand(n,1));
 x = x0 + r.*cos(t);
 y = y0 + r.*sin(t);
 

n1=5*n;R1=0.75;R2=1/R1;x01=1; y01=1.5+R2-1;
rng(n1)
t = 2*pi*rand(n1,1);
r1 = R1*sqrt(rand(n1,1)); r2 = r1*R2/R1;
x1 = x01 + r1.*cos(t);
y1 = y01 + r2.*sin(t);
xy1_rec=[]; 
 for i=1:n1
 if ((x1(i)-1)^2+(y1(i)-1)^2)^0.5>R
      xy1_rec=[xy1_rec;x1(i),y1(i)];
 end
 end
n2=3*n;
location=[x,y];
location2=xy1_rec(1:3*n,:);

connection_d=0;
for i=1:n
for j=1:n2
    if ((location2(j,1)-location(i,1))^2+(location2(j,2)-location(i,2))^2)^0.5<dmax
        connection_d=connection_d+1;
    end
end
end
connection_d_ave=connection_d/n2;
end
%% p
function connection_p_ave=cal_connection_p(p)
n=60; 
R=0.5; x0=1; y0=1;
rng(n)
 t = 2*pi*rand(n,1);
 r = R*sqrt(rand(n,1));
 x = x0 + r.*cos(t);
 y = y0 + r.*sin(t);
 

n1=5*n;R1=0.75;R2=1/R1;x01=1; y01=1.5+R2-1;
rng(n1)
t = 2*pi*rand(n1,1);
r1 = R1*sqrt(rand(n1,1)); r2 = r1*R2/R1;
x1 = x01 + r1.*cos(t);
y1 = y01 + r2.*sin(t);
xy1_rec=[]; 
 for i=1:n1
 if ((x1(i)-1)^2+(y1(i)-1)^2)^0.5>R
      xy1_rec=[xy1_rec;x1(i),y1(i)];
 end
 end
n2=3*n;
location=[x,y];
location2=xy1_rec(1:3*n,:);

connection_p=0;
for i=1:n
for j=1:n2
    xx=rand;
if xx<p
  connection_add=1;
else
  connection_add=0;
end
        connection_p=connection_p+connection_add;
   
end
end
connection_p_ave=connection_p/n2;
end
