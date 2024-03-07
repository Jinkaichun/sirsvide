function [obj_gen,relation_matrix,pa_matrix]=getrelation(strains,pa_strains,gen_group)
n_s=0;
for i=(1:length(gen_group))
    gen=gen_group{i};
    gen=gen(:);
    gen(find(gen==0))=[];
    if length(gen)>n_s
        max_gen=i;
        n_s=length(gen);
    end
end
obj_gen=gen_group{max_gen};
obj_gen=obj_gen(:);
obj_gen(find(obj_gen==0))=[];
n_gen=length(obj_gen);
relation_matrix=zeros(n_gen,n_gen);
pa_matrix=zeros(n_gen,n_gen);
for i=(1:n_gen)
    for j=(1:n_gen)
        s_a=obj_gen(i);
        s_b=obj_gen(j);
        p_a=s_a;
        p_b=s_b;
        relation=0;
        while p_a~=p_b
            p_a=pa_strains(find(strains==p_a));
            p_b=pa_strains(find(strains==p_b));
            relation=relation+1;
        end
        relation_matrix(i,j)=relation;
        pa_matrix(i,j)=p_a;
    end
end
end