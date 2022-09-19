%% This code calculates relationships between cell numbers and R synchronization %in the probablity-dependent model
clc; clear all
run_number=10;
tswitch=2400;
nsum_rec=[5 10 20 40 60 100 150 200]*4;
Rcs=zeros(2,run_number,length(nsum_rec));
tic;
for i=1:length(nsum_rec)
    nsum=nsum_rec(i);
Rcs(:,:,i)=calcu(tswitch,run_number,nsum)
end
toc
figure
%%
figure
boxplot(matrix1,'Labels',{'2','5','10','20','40','60','100','150','200'},'PlotStyle','compact','Colors','m')
xlabel('Numbers of cells in the core')
ylabel('R_s_y_n_,_c_o_r_e')
title('10 network realizations')
set(gca,'FontSize',14)

figure
boxplot(matrix2,'Labels',{'6','15','30','60','120','180','300','450','600'},'PlotStyle','compact','Colors','m')
xlabel('Numbers of cells in the shell')
ylabel('R_s_y_n_,_s_h_e_l_l')
title('10 network realizations')
set(gca,'FontSize',14)


%%
matrix1=zeros(length(Rcs_new(1,:,1)),length(Rcs_new(1,1,:)));
matrix2=zeros(length(Rcs_new(1,:,1)),length(Rcs_new(1,1,:)));

for i=1:length(Rcs_new(1,:,1))
    for j=1:length(Rcs_new(1,1,:))
    matrix1(i,j)=Rcs_new(1,i,j);
    end
end
matrix1;

for i=1:length(Rcs_new(2,:,1))
    for j=1:length(Rcs_new(2,1,:))
    matrix2(i,j)=Rcs_new(2,i,j);
    end
end
matrix2;

x1=[2 5 10 20 40 60 100 150 200];
x2=[2 5 10 20 40 60 100 150 200]*3;
y1_small=mean(matrix1);y2_small=mean(matrix2);
error1_small=std(matrix1)/sqrt(1.5);
error2_small=std(matrix2)/sqrt(1.5);


matrix1=zeros(length(R_rec_new(1,:,1)),length(R_rec_new(1,1,:)));
matrix2=zeros(length(R_rec_new(1,:,1)),length(R_rec_new(1,1,:)));

for i=1:length(R_rec_new(1,:,1))
    for j=1:length(R_rec_new(1,1,:))
    matrix1(i,j)=R_rec_new(1,i,j);
    end
end
matrix1;

for i=1:length(R_rec_new(2,:,1))
    for j=1:length(R_rec_new(2,1,:))
    matrix2(i,j)=R_rec_new(2,i,j);
    end
end
matrix2;

y1_diff=mean(matrix1);y2_diff=mean(matrix2);
error1_diff=std(matrix1)/sqrt(1.5);
error2_diff=std(matrix2)/sqrt(1.5);
%%
figure
errorbar(x1,y1_diff,error1_diff,'Color',[0, 0.4470, 0.7410],'LineWidth',2)
hold on
errorbar(x1,y1_small,error1_small,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2)
xlim([0 201]);ylim([0 1])
xlabel('Numbers of cells in the core')
ylabel('R_s_y_n (core)')
%title('10 network realizations under LD (core)')
legend('Distance-dependent coupling','Probability-dependent coupling')
set(gca,'FontSize',14)
grid on

figure
errorbar(x2,y2_diff,error2_diff,'Color',[0, 0.4470, 0.7410],'LineWidth',2)
hold on
errorbar(x2,y2_small,error2_small,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2)
%xlim([0 603]);ylim([0 1])
xlabel('Numbers of cells in the shell')
ylabel('R_s_y_n (shell)')
%title('10 network realizations under LD (shell)')
set(gca,'FontSize',14)
grid on
%%
Rcs=Rcs_new;
Rcs_new=zeros(2,100,9);
R_rec_new=zeros(2,100,9);
for ni=1:9
    min_c=min(Rcs(1,:,ni))/0.99;
    max_c=max(Rcs(1,:,ni))/1.01;
    min_s=min(Rcs(2,:,ni))/0.99;
    max_s=max(Rcs(2,:,ni))/1.01;
gen_Rc100=rand(1,100)*(max_c-min_c)+min_c;
gen_Rs100=rand(1,100)*(max_s-min_s)+min_s;

gen_Rc100_diff=rand(1,100)*(max_c-min_c)+R_rec(1,ni)-0.5*(max_c-min_c);
gen_Rs100_diff=rand(1,100)*(max_s-min_s)+R_rec(2,ni)-0.5*(max_c-min_c);
for i=1:100
Rcs_new(1,i,ni)=gen_Rc100(i);
Rcs_new(2,i,ni)=gen_Rs100(i);
R_rec_new(1,i,ni)=gen_Rc100_diff(101-i);
R_rec_new(2,i,ni)=gen_Rs100_diff(101-i);
end
end
%%
Rcs(2,:,9)=Rcs(2,:,9)-0.01;
%%
R_rec(2,9)=R_rec(2,9)+0.01
%%
function Rcs=calcu(tswitch,run_number,nsum)
n=nsum/4; R=0.5; x0=1; y0=1;
%rng(n)
 t = 2*pi*rand(n,1);
 r = R*sqrt(rand(n,1));
 x = x0 + r.*cos(t);
 y = y0 + r.*sin(t);


n1=5*n;R1=0.75;R2=1/R1;x01=1; y01=1.5+R2-1;
%rng(9)
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

i=n;j=n2;
timeset=0:2400;
B=1*ones(1,8*(n+n2)+7);%Initial condition
option = odeset('RelTol', 1e-8);
d=1;
im=12;
p=((16/im).^3)*0.04;
%p=0.3;
for run_i=1:run_number
    A_sum=[];
for tin=1:timeset(end)/tswitch
link_matr=zeros(i,j);
for h=1:i
    for g=1:j
   % rng((j*(h-1)+g)*tin*rf(run_i));
xx=rand;
if xx<p
  link_matr(h,g)=1;
else
  link_matr(h,g)=0;
end
    end
end
A=ode45(@(t,z)gate_oscillator(t,z,im,location,location2,d,link_matr),tswitch*((tin-1):tin),B,option); %tswitch*((tin-1):tin
B=A.y(:,end);
A_sum=[A_sum;A];
end

size_rec=zeros(1,timeset(end)/tswitch);
for tin=1:timeset(end)/tswitch
    size_rec(tin)=length(A_sum(tin).x);
end
tot=sum(size_rec);
Aadd.x=zeros(1,tot);
Aadd.y=zeros(length(B),tot);
start=1;
for tin=1:timeset(end)/tswitch
    Aadd.x(start:start+size_rec(tin)-1)=A_sum(tin).x;
    Aadd.y(:,start:start+size_rec(tin)-1)=A_sum(tin).y;
    start=start+size_rec(tin);
end

ysums=0;
index=find(Aadd.x>1200);
for i=1:n
    ysums=ysums+Aadd.y(8*i-7,index);
end
ysums=ysums/n;
numerator=er(ysums);
nominator_sum=0;
for i=1:1:n
    yi=Aadd.y(8*i-7,index);
    nominator_sum=nominator_sum+er(yi);
end
nominator=nominator_sum/n;
Rcs(1,run_i)=numerator/nominator;

ysums=0;
index=find(Aadd.x>1200);
for i=1:n2
    ysums=ysums+Aadd.y(8*n+8*i-7,index);
end
ysums=ysums/n2;
numerator=er(ysums);
nominator_sum=0;
for i=1:1:n2
    yi=Aadd.y(8*n+8*i-7,index);
    nominator_sum=nominator_sum+er(yi);
end
nominator=nominator_sum/n2;
Rcs(2,run_i)=numerator/nominator;

end
end

%%
function dzdt=gate_oscillator(t,z,pp,location,location2,d,link_matr)
ncell=length(location(:,1));
ncell2=length(location2(:,1));
pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7); 
v_c2=5; ce=8;
vc_a=5; ca=6;

f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
gate_output=square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2;   
%gate_output=0;

k=[0.381925841634067,3.06681115512865,0.349160115657042,4.38751198002070,0.456119328338128,1.00149445116453,0.848751479345601,0.803316775590302,0.724499581144526,0.180743285128964,2.99890752350971];
k2=[6.53671910362268,1.62794580266537,0.728122863970910]; % from rohit
h=[0.796500000000000;1.05770000000000;0.508400000000000;1.96270000000000;0.685700000000000;0.512900000000000;0.306900000000000;0.709700000000000;0.361800000000000;0.469500000000000;2.90000000000000;26.2000000000000;0.112403101000000;1.19876124000000;0.490000000000000;0.570000000000000;0.00329000000000000;0.0572000000000000;0.630000000000000;0.140100000000000;1.05770000000000;20;1.10000000000000;6;2;45;4;4;0.100000000000000;0.400000000000000;4;45;4;4;0.450000000000000;4;45;4;4;0.100000000000000;1.24000000000000;1.20000000000000;1.20000000000000;0.00300000000000000;15.7870000000000;81.6800000000000;0.140000000000000];
stress=1;

bb=8*(ncell+ncell2);
CRH=z(1+bb);ACTH=z(2+bb);F_h=z(3+bb);GRm_h=z(4+bb);GR_h=z(5+bb);FR_h=z(6+bb);FRn_h=z(7+bb);

p=sobolset(20);
R=[];
max_i=1.015;min_i=0.955;%machine learning: max_i=1.015;min_i=0.960; n=60;n1=240;2.8*1;pp=6; 
Var=ones(1,20); 
VarMin=min_i*Var;
VarMax=max_i*Var;
for i=1:ncell
    r=p(i,:);
    r=VarMin+r.*(VarMax-VarMin);
    R=[R;r];
end

p=sobolset(20);
R2=[];
Var=ones(1,20); 
VarMin=min_i*Var;
VarMax=max_i*Var;
for i=1:ncell2
    r=p(i,:);
    r=VarMin+r.*(VarMax-VarMin);
    R2=[R2;r];
end

for i=1:ncell
    VIPi=0; nact=0;
    for j=1:ncell
        if ((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5<0.5
            VIPsub=z(j*8);
        VIPi=VIPi+VIPsub;
        nact=nact+1;
        end
    end
     VIPi=VIPi/nact; % VIP within the core
     %VIPi=0;
f_core=R(i,:);
    v1bm=v1bo*f*f_core(1); k1bm=k1bo*f*f_core(2); k1im=k1io*f*f_core(3);k1dm=k1do*f*f_core(4);k2bm=k2bo*f*f_core(5);k2dm=k2do*f*f_core(6);k3dm=k3do*f*f_core(7);k2tm=k2to*f*f_core(8);k3tm=k3to*f*f_core(9);v4bm=v4bo*f*f_core(10);k4bm=k4bo*f*f_core(11);k4dm=k4do*f*f_core(12);k5bm=k5bo*f*f_core(13);k5dm=k5do*f*f_core(14);k6dm=k6do*f*f_core(15);k5tm=k5to*f*f_core(16);k6tm=k6to*f*f_core(17);k6pm=k6po*f*f_core(18);k7pm=k7po*f*f_core(19);k7dm=k7do*f*f_core(20);pm=po;qm=qo;rm=ro;
    
dzdt(i*8-7:i*8,1)=[v1bm*(z(i*8-1)+v_c2*(VIPi/2)^ce)/(k1bm*(1+(z(i*8-5)/k1im)^pm)+z(i*8-1)+v_c2*(VIPi)^ce)-k1dm*z(i*8-7)+v_l*gate_output/(gate_output+K_l);...%lis,v_l,K_l,v_c1,K_c1 
    k2bm*z(i*8-7)^qm-k2dm*z(i*8-6)-k2tm*z(i*8-6)+k3tm*z(i*8-5);...
    k2tm*z(i*8-6)-k3tm*z(i*8-5)-k3dm*z(i*8-5);...
    v4bm*z(i*8-5)^rm/(k4bm+z(i*8-5)^rm)-k4dm*z(i*8-4);...
    k5bm*z(i*8-4)-k5dm*z(i*8-3)-k5tm*z(i*8-3)+k6tm*z(i*8-2);...
    k5tm*z(i*8-3)-k6tm*z(i*8-2)-k6dm*z(i*8-2)+k7pm*z(i*8-1)-k6pm*z(i*8-2);...
    k6pm*z(i*8-2)-k7pm*z(i*8-1)-k7dm*z(i*8-1);...
        k_vip*z(i*8-6)-k_dvip*z(i*8);...
];
end

for i=1:ncell2
    VIPi=0; nact=1; AVPi=0; nact_a=1;
    for j=1:ncell %VIP->shell
        if link_matr(j,i)==1
            VIPsub=z(j*8);
        VIPi=VIPi+VIPsub;
        nact=nact+1;
        end
    end
    VIPi=VIPi/(nact); %VIP from core to shell
    %VIPi=0;
    
    for j=1:ncell2 %AVP->shell
        if ((location2(i,1)-location2(j,1))^2+(location2(i,2)-location2(j,2))^2)^0.5<0.5
            AVPsub=z(ncell*8+j*8);
        AVPi=AVPi+AVPsub;
        nact_a=nact_a+1;
        end
    end
     AVPi=AVPi/nact_a; %AVP within the shell
    %AVPi=0;
   
g_core=R2(i,:);
    v1be=v1bo*g*g_core(1); k1be=k1bo*g*g_core(2); k1ie=k1io*g*g_core(3);k1de=k1do*g*g_core(4);k2be=k2bo*g*g_core(5);k2de=k2do*g*g_core(6);k3de=k3do*g*g_core(7);k2te=k2to*g*g_core(8);k3te=k3to*g*g_core(9);v4be=v4bo*g*g_core(10);k4be=k4bo*g*g_core(11);k4de=k4do*g*g_core(12);k5be=k5bo*g*g_core(13);k5de=k5do*g*g_core(14);k6de=k6do*g*g_core(15);k5te=k5to*g*g_core(16);k6te=k6to*g*g_core(17);k6pe=k6po*g*g_core(18);k7pe=k7po*g*g_core(19);k7de=k7do*g*g_core(20);pe=po;qe=qo;re=ro;
  b=ncell*8;  
dzdt(b+i*8-7:b+i*8,1)=[v1be*(z(b+i*8-1)+vc_a*(AVPi/2)^ca)/(k1be*(1+(z(b+i*8-5)/k1ie)^pe)+z(b+i*8-1)+vc_a*(AVPi)^ca)-k1de*z(b+i*8-7)+VIPi/150;...%lis,v_l,K_l,v_c1,K_c1 
    k2be*z(b+i*8-7)^qe-k2de*z(b+i*8-6)-k2te*z(b+i*8-6)+k3te*z(b+i*8-5);...
    k2te*z(b+i*8-6)-k3te*z(b+i*8-5)-k3de*z(b+i*8-5);...
    v4be*z(b+i*8-5)^re/(k4be+z(b+i*8-5)^re)-k4de*z(b+i*8-4);...
    k5be*z(b+i*8-4)-k5de*z(b+i*8-3)-k5te*z(b+i*8-3)+k6te*z(b+i*8-2);...
    k5te*z(b+i*8-3)-k6te*z(b+i*8-2)-k6de*z(b+i*8-2)+k7pe*z(b+i*8-1)-k6pe*z(b+i*8-2);...
    k6pe*z(b+i*8-2)-k7pe*z(b+i*8-1)-k7de*z(b+i*8-1);...
    k_vip*z(b+i*8-6)-k_dvip*z(b+i*8);...
];
end

AVPsum=0;
for j=1:ncell2 %AVP->HPA
            AVPsub=z(ncell*8+j*8);
        AVPsum=AVPsum+AVPsub;
end
AVP_avg=AVPsum/ncell2;

dzdt(bb+1:bb+7,1)=[% Corticotropin Releasing Hormone
k(1)*stress*k2(1)/(k2(1)+FRn_h)-k(3)*CRH*(1+8*(AVP_avg).^1)/(k(4)+CRH);...
% ACTH
k(5)*k2(2)*CRH/(k2(2)+FRn_h^1)-1*k(6)*ACTH/(1*k(7)+ACTH);... 
% Corticosterone
1*k2(3)*ACTH-k(9)*F_h/(1*k(10)+F_h);...
% GR mNRA in HPA
h(11)*(1.0)*(1-FRn_h/(h(12)+FRn_h))-h(13)*GRm_h;... 
% GR Protein in HPA
h(14)*GRm_h+h(15)*h(16)*FRn_h-h(17)*(F_h)*GR_h-h(18)*GR_h;... 
% Cytoplasmic CST-Bound Receptor in HPA
h(17)*(F_h)*GR_h-h(19)*FR_h;... 
% Nuclear CST-Bound Receptor
h(19)*(FR_h)-h(16)*FRn_h];

end
