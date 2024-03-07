clc,clear,close all
%-----parameters set-----
N=10^9;
gam=1/10;
R0_0=3;
patho_0=0.5;%initial pathogenicity
b=R0_0*gam;
T=90;%time interval between reinfection
omega=1/T;
Ti=500;%time of humoral immunity waning
omegai=1/Ti;
m=10^-8; %mutation rate per infected per day
aa=300;
genome_0=zeros(1,aa);%initial genome
t_total=5000;%simulation time(day)
dt=1;
dev_R0=1;
dev_patho=0.1;
n_mutsite=3;
parameter_set=[N,gam,R0_0,patho_0,T,Ti,m,aa,t_total,dt,dev_R0,dev_patho,n_mutsite];
total_rep=1;
escape_score=ones(1,100);
load('escape_score')
fit=0.1;%probability of beneficial mutation
k=0.5;%contact rate attenuation

%-----recording variables-----
It_r=cell(total_rep,1);
St_r=cell(total_rep,1);
Dt_r=cell(total_rep,1);
Rt_r=cell(total_rep,1);
Sit_r=cell(total_rep,1);
casest_r=cell(total_rep,1);
strains_remaint_r=cell(total_rep,1);
pa_strains_r=cell(total_rep,1);
R0_strains_r=cell(total_rep,1);
T0_strains_r=cell(total_rep,1);
patho_strains_r=cell(total_rep,1);
R0_ave_r=cell(total_rep,1);
patho_ave_r=cell(total_rep,1);
shannon_r=cell(total_rep,1);
pop_genome_r=cell(total_rep,1);
pop_escape_r=cell(total_rep,1);

%-----replicate------
tic
for rep=(1:total_rep)
    rep
    %-----initial value-----
    pop_genome=genome_0;%genome of population
    I=10;
    S=N-I;
    R=0;
    Si=0;
    D=0;
    imatrix=0;
    R0_strains=R0_0;%R_0 of strains
    patho_strains=patho_0;
    pa_strains=0;%parent of strains
    T0_strains=0;%birth time of strains
    b_strains=R0_strains*gam;
    strains_remain=1;
    strain_max=100;
    t=0;
    num=0;
    It=[];
    Sit=[];
    Dt=[];
    Rt=[];
    St=[];
    casest=[];
    strains_remaint=[];
    R0_ave=[];
    patho_ave=[];
    shannon=[];
    pop_escape=[];

    %-----simulation------
    while t<t_total
        num=num+1;
        %-----get dX-----
        imatrix_remain=imatrix(strains_remain,strains_remain);
        dI=exp(-k*patho_strains(strains_remain)).*b_strains(strains_remain).*I*S/N-gam*I+(imatrix_remain*Si).*exp(-k*patho_strains(strains_remain)).*b_strains(strains_remain).*I/N;
        dS=-sum(exp(-k*patho_strains(strains_remain)).*b_strains(strains_remain).*I*S/N)+sum(omegai*Si);
        dR=gam*(1-0.1*patho_strains(strains_remain)).*I-omega*R;
        dSi=omega*R-omegai*Si-imatrix_remain*(I.*exp(-k*patho_strains(strains_remain)).*b_strains(strains_remain)).*Si/N;
        dD=sum(gam*0.1*patho_strains(strains_remain).*I);
        cases=exp(-k*patho_strains(strains_remain)).*b_strains(strains_remain).*I*S/N+(imatrix_remain*Si).*exp(-k*patho_strains(strains_remain)).*b_strains(strains_remain).*I/N;
        %-----update-----
        I=max(I+dI*dt,0);
        S=max(S+dS*dt,0);
        R=max(R+dR*dt,0);
        Si=max(Si+dSi*dt,0);
        D=max(D+dD*dt,0);
        I=I.*(I>=1);
        del_index=find(I<1&R<1000&Si<1000);
        rem_index=find(I>=1|R>=1000|Si>=1000);
        strains_remain=strains_remain(rem_index);
        cases(del_index)=[];
        I(del_index)=[];
        Si(del_index)=[];
        R(del_index)=[];
        I_sum=sum(I);
        R_sum=sum(R);
        
        %-----mutation-----
        n_strains=length(R0_strains);
        nmut=random('Poisson',I*m);
        for i=(1:length(I))
            if nmut(i)>0
                mutnum=random('Poisson',n_mutsite*ones(nmut(i),1));
                %mutnum=n_mutsite*ones(nmut(i),1);
                genome=pop_genome(strains_remain(i),:);
                size_200=[];
                genome_new=[];
                for j=(1:nmut(i))
                    seed=randperm(aa);
                    mutsite=seed(1:mutnum(j));
                    size_200(j,:)=length(find(mutsite<201));
                    newaa=randi([0,19],1,mutnum(j));
                    geo=genome;
                    geo(mutsite)=newaa;
                    genome_new(j,:)=geo;
                end
                I(i)=I(i)-nmut(i);
                %R0_new=R0_strains(strains_remain(i))+size_200.*(gamrnd(ones(nmut(i),1),R0_strains(strains_remain(i)))-R0_strains(strains_remain(i)));   %根据gamma分布生成变异株R0
                seed=rand(nmut(i),1);
                fitting=zeros(nmut(i),1);
                ind_fit=find(seed<fit);
                ind_del=find(seed>=fit);
                fitting(ind_fit)=1;
                fitting(ind_del)=-1;
                R0_new=R0_strains(strains_remain(i))+fitting.*size_200.*abs(randn(nmut(i),1))*dev_R0;
                patho_new=patho_strains(strains_remain(i))+size_200.*randn(nmut(i),1)*dev_patho;%根据正态分布生成变异株patho
                R0_new=(R0_new>0).*R0_new;%R_0 is positive
                patho_new=(patho_new>0).*patho_new;
                %patho_new(find(patho_new>1))=1;
                R0_strains=cat(1,R0_strains,R0_new);
                patho_strains=cat(1,patho_strains,patho_new);
                pop_genome=cat(1,pop_genome,genome_new);
                pa_strains=cat(1,pa_strains,strains_remain(i)*ones(nmut(i),1)); 
                T0_strains=cat(1,T0_strains,t*ones(nmut(i),1)); 
                I=cat(1,I,ones(nmut(i),1));
                Si=cat(1,Si,zeros(nmut(i),1));
                R=cat(1,R,zeros(nmut(i),1));
                cases=cat(1,cases,ones(nmut(i),1));
                strains_remain=cat(1,strains_remain,(length(R0_strains)-nmut(i)+1:length(R0_strains))');
            end
        end
        b_strains=R0_strains*gam;
        %-----update immune matrix-----
        n_strains_new=length(R0_strains)-n_strains;
        if n_strains_new>0
            imatrix_new=zeros(length(R0_strains),length(R0_strains));
            imatrix_new(1:n_strains,1:n_strains)=imatrix;
            plus=zeros(length(R0_strains),n_strains_new);
            for i=(1:n_strains_new)
                plus(:,i)=align(n_strains+i,pop_genome(:,1:100),escape_score);
            end
            imatrix_new(:,n_strains+1:end)=plus;
            imatrix_new(n_strains+1:end,:)=plus';
            imatrix=imatrix_new;
        end
        %-----record-----
        It(1:length(I),num)=I;
        casest(1:length(I),num)=cases;
        strains_remaint(1:length(I),num)=strains_remain;
        St(num)=S;
        Dt(num)=D;
        Rt(1:length(R),num)=R;
        Sit(1:length(Si),num)=Si;
        I_sumt(num)=I_sum;
        R_sumt(num)=R_sum;
        R0_ave(num)=sum(I.*R0_strains(strains_remain)/I_sum);
        patho_ave(num)=sum(I.*patho_strains(strains_remain)/I_sum);
        P=I/I_sum;
        P(find(P==0))=[];
        shannon(num)=-sum(P.*log2(P));
        imatrix_remain=imatrix(strains_remain,strains_remain);
        Si_p=Si/sum(Si);
        I_p=I/sum(I);
        strains_escape=sum(Si_p.*imatrix_remain);
        pop_escape(num)=sum(strains_escape.*I_p');
        t=t+dt;
    end
    It_r{rep}=It;
    Rt_r{rep}=Rt;
    St_r{rep}=St;
    Dt_r{rep}=Dt;
    Sit_r{rep}=Sit;
    casest_r{rep}=casest;
    strains_remaint_r{rep}=strains_remaint;
    pa_strains_r{rep}=pa_strains;
    R0_strains_r{rep}=R0_strains;
    T0_strains_r{rep}=T0_strains;
    patho_strains_r{rep}=patho_strains;
    R0_ave_r{rep}=R0_ave;
    patho_ave_r{rep}=patho_ave;
    shannon_r{rep}=shannon;
    pop_genome_r{rep}=pop_genome;
    pop_escape_r{rep}=pop_escape;
end
toc
%%
filename='test_result.mat';

save(filename,'parameter_set','escape_score','It_r','St_r','Dt_r','Rt_r','Sit_r','casest_r','strains_remaint_r','pa_strains_r','R0_strains_r','T0_strains_r','patho_strains_r','R0_ave_r','patho_ave_r','shannon_r','pop_genome_r','pop_escape_r')

