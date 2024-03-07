function [conver_matrix,conver_freq]=getconvergent(pop_genome,pa_matrix)
[a,b]=size(pa_matrix);
conver_matrix=cell(a,a);
conver_sum=zeros(1,300);
for i=(1:a)
    for j=(1:b)
        s_a=pop_genome(pa_matrix(i,i),:);
        s_b=pop_genome(pa_matrix(j,j),:);
        s_p=pop_genome(pa_matrix(i,j),:);
        conver=(s_a~=0)&(s_b~=0)&(s_p==0);
        conver_matrix{i,j}=conver;
        conver_sum=conver_sum+conver;
    end       
end
conver_freq=conver_sum/(a^2-a);
end