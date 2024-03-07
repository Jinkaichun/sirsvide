clc,clear,close all

filename='test_result.mat';
load(filename)
%%
close all
%-----data process-----

%-----load data-----
rep=1;
dicname=['rep_',num2str(rep)];

It=It_r{rep};
Rt=Rt_r{rep};
St=St_r{rep};
Dt=Dt_r{rep};
Sit=Sit_r{rep};
casest=casest_r{rep};
strains_remaint=strains_remaint_r{rep};
pa_strains=pa_strains_r{rep};
R0_strains=R0_strains_r{rep};
T0_strains=T0_strains_r{rep};
patho_strains=patho_strains_r{rep};
R0_ave=R0_ave_r{rep};
patho_ave=patho_ave_r{rep};
shannon=shannon_r{rep};
pop_genome=pop_genome_r{rep};
pop_escape=pop_escape_r{rep};

%-----load parameters------
N=parameter_set(1);
gam=parameter_set(2);
R0_0=parameter_set(3);
patho_0=parameter_set(4);
T=parameter_set(5);
Ti=parameter_set(6);
m=parameter_set(7);
aa=parameter_set(8);
t_total=parameter_set(9);
dt=parameter_set(10);
tt=(dt:dt:t_total);

I_limit=N*0.001;%定义VOC的阈值
%-----get sum-----
n_strains=length(R0_strains);
if n_strains>1
    I_sumt=sum(It);
    R_sumt=sum(Rt);
    Si_sumt=sum(Sit);
else
    I_sumt=It;
    R_sumt=Rt;
    Si_sumt=Sit;
end
%-----get imatrix-----

imatrix=zeros(n_strains,n_strains);
for i=(1:n_strains)
    imatrix(i,:)=align(i,pop_genome(:,1:100),escape_score);
end

It_all=[];
Rt_all=[];
Sit_all=[];
%-----remain 2 all------
for i=(1:length(tt))
    It_all_i=zeros(length(R0_strains),1);
    Rt_all_i=zeros(length(R0_strains),1);
    Sit_all_i=zeros(length(R0_strains),1);
    strains_remain=strains_remaint(:,i);
    n_index=sum(strains_remain>0);
    strains_remain(n_index+1:end)=[];
    It_all_i(strains_remain)=It(1:n_index,i);
    Rt_all_i(strains_remain)=Rt(1:n_index,i);
    Sit_all_i(strains_remain)=Sit(1:n_index,i);
    It_all(:,i)=It_all_i;
    Rt_all(:,i)=Rt_all_i;
    Sit_all(:,i)=Sit_all_i;
end

%-----get dominant strains-----

if n_strains>1
    I_max=max(It_all,[],2);
    del_index=find(I_max<I_limit);
    strains_main=find(I_max>=I_limit);
    It_main=It_all(strains_main,:);
    It_others=I_sumt-sum(It_main,1);
    peak=max(I_sumt);
    [domi_num,domi_strain]=max(It_all);
else
    I_max=max(It_all);
    strains_main=1;
    peak=max(I_sumt);
    It_main=It_all;
    It_others=zeros(1,t_total);
    domi_strain=ones(1,t_total);
end

domi_R0=R0_strains(domi_strain);
%domi_escape=R0_strains(domi_strain);
domi_patho=patho_strains(domi_strain);


%-----infection statement-----


colormap hsv
cmap = colormap;
figure(1)
subplot(2,1,1)
pp=area(tt,[It_main;It_others]');
n_It_main=size(It_main,1);
load('color_order.mat')
color_order=cmap(order(1:n_It_main),:);
for i=(1:n_It_main)
    pp(i).FaceColor=color_order(i,:);
end

pp(end).FaceColor=[0.5,0.5,0.5];

axis([0,t_total,0,1.1*peak])
xlabel('Time(days)','FontWeight','bold')
ylabel('Number of infections','FontWeight','bold')


%-----trait evolution-----
figure(2)

hold on
f_alpha=0.4;
scatter(T0_strains,R0_strains*gam,'filled','b','MarkerFaceAlpha',f_alpha)
plot(tt,R0_ave*gam,'LineWidth',2,'Color','r')
xlabel('Time(days)','FontWeight','bold')
ylabel('Transmissibility(\beta)','FontWeight','bold')

legend('\beta of every strains','Average \beta')


figure(3)
hold on
scatter(T0_strains,patho_strains,'filled','b','MarkerFaceAlpha',f_alpha)
plot(tt,patho_ave,'LineWidth',2,'Color','r')
xlabel('Time(days)','FontWeight','bold')
ylabel('Pathogenicity(\alpha)','FontWeight','bold')
ylim([0,1.1])
legend('\alpha of every strains','Average \alpha')

figure(4)
plot(tt,pop_escape,'LineWidth',2,'Color','r')
xlabel('Time(days)','FontWeight','bold')
ylabel('Population immune escape capability','FontWeight','bold')

%-----alpha diversity-----
figure(5)
plot(tt,shannon,'LineWidth',2,'Color','r')
xlabel('Time(days)','FontWeight','bold')
ylabel('Strains diversity (Shannon index)','FontWeight','bold')



%-----draw tree-----
pa_strains_main=pa_strains(strains_main);
C=setdiff(pa_strains_main,strains_main);
C(1)=[];
while ~isempty(C)
    strains_main=cat(1,strains_main,C);
    pa_strains_main=pa_strains(strains_main);
    C=setdiff(pa_strains_main,strains_main);
    C(1)=[];
end
n_It_main=length(strains_main);
T0_strains_main=T0_strains(strains_main);
I_max_main=I_max(strains_main);
R0_strains_main=R0_strains(strains_main);
patho_strains_main=patho_strains(strains_main);
color_order2=cmap(order(1:n_It_main),:);

for i=(1:n_It_main)
    T_ap=T0_strains_main(i);
    num=find(tt==T_ap);
    Si=Sit_all(:,num);
    Sip=Si/sum(Si);
    diversity_i=imatrix(:,strains_main(i));
    IE_strains_m(i)=sum(diversity_i.*Sip);
end
IE_strains_main=IE_strains_m';
figure(6)
[gen_group,X_all,Y_all]=drawtree_gen(strains_main,pa_strains_main,T0_strains_main,I_max_main,color_order2,R0_strains_main*gam,patho_strains_main,IE_strains_main);
[obj_gen,relation_matrix,pa_matrix]=getrelation(strains_main,pa_strains_main,gen_group);
[conver_matrix,conver_freq]=getconvergent(pop_genome,pa_matrix);

%%
n_VOC=length(strains_main)-1;
VOC=strains_main(2:end);
VOC_pa=pa_strains_main(2:end);
R0_VOC=R0_strains(VOC);
patho_VOC=patho_strains(VOC);
R0_VOC_pa=R0_strains(VOC_pa);
patho_VOC_pa=patho_strains(VOC_pa);
R0_evo_trend=R0_VOC-R0_VOC_pa;
patho_evo_trend=patho_VOC-patho_VOC_pa;





function ave=getave(X,T)
ave=reshape(X,[T,length(X)/T]);
ave=mean(ave);
end