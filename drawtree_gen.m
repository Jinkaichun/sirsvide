function [gen_group,X_all,Y_all]=drawtree_gen(strains,parents,T0,I_max,color_order,R0,patho,IE) %输入等长向量
%----generation group------
gen=0;
node=strains(find(parents==0));
gen_group={};
gen_group{gen+1}=node;
sons=node;
n_sons=length(sons);
Y_all=linspace(0,100,n_sons+2);
Y_all=Y_all(2:end-1);
X=T0(find(parents==0));
order=node;
gen_nodes(gen+1)=n_sons;
size_all=I_max(find(parents==0));
color_all=color_order(find(parents==0),:);
R0_all=R0(find(parents==0));
patho_all=patho(find(parents==0));
IE_all=IE(find(parents==0));
while ~isempty(sons)
    gen=gen+1;
    node=[];
    last_gen=gen_group{gen};
    last_gen=last_gen(:);
    last_gen(find(last_gen==0))=[];
    num=0;
    for pa=last_gen'
        num=num+1;
        son_index=find(parents==pa);
        son=strains(son_index);%找到上一世代每个节点的子节点
        if isempty(son)
            node(1,num)=0;
        else
            son_T0=T0(son_index);
            son_size=I_max(son_index);
            son_color=color_order(son_index,:);
            son_R0=R0(son_index);
            son_patho=patho(son_index);
            son_IE=IE(son_index);
            [T,I]=sort(son_T0,'descend');
            node(1:length(son),num)=son(I);
            size=son_size(I);
            size_all=cat(1,size_all,size);
            color=son_color(I,:);
            color_all=cat(1,color_all,color);
            R00=son_R0(I);
            R0_all=cat(1,R0_all,R00);
            Pa=son_patho(I);
            patho_all=cat(1,patho_all,Pa);
            ie=son_IE(I);
            IE_all=cat(1,IE_all,ie);
%             X=cat(1,X,T);
        end
    end
    gen_group{gen+1}=node;
    sons=node(:);
    sons(find(sons==0))=[];
    n_sons=length(sons);
    if ~isempty(sons)
        order=cat(1,order,sons);
        Y=linspace(100,0,n_sons+2);
        Y=Y(2:end-1)';
        Y_all=cat(1,Y_all,Y);
    end
    gen_nodes(gen+1)=n_sons;
end
X=linspace(0,100,gen+3);
X=X(2:end-1)';
X_all=[];
for i=(1:length(gen_group)-1)
    X_all=cat(1,X_all,X(i)*ones(gen_nodes(i),1));
end
bubblechart(X_all,Y_all,size_all,color_all)
axis([0,100,0,100])
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
hold on
R0_str={};
patho_str={};
IE_str=[];
for i=(1:length(R0_all))
    R0_str{i}=num2str(round(R0_all(i),2));
    patho_str{i}=num2str(round(patho_all(i),2));
    IE_str{i}=num2str(round(IE_all(i),2));
end
fsize=7;
text(X_all,Y_all+1.6,R0_str,'HorizontalAlignment','center','FontSize',fsize,'FontWeight','bold')
text(X_all,Y_all-1.5,patho_str,'HorizontalAlignment','center','FontSize',fsize,'FontWeight','bold')
% text(X_all,Y_all-2.2,IE_str,'HorizontalAlignment','center','FontSize',fsize)
for i=(1:length(gen_group)-2)
    gen=gen_group{i};
    gen=gen(:);
    gen(find(gen==0))=[];
    next_gen=gen_group{i+1};
    for j=(1:length(gen))
        index_p=find(order==gen(j));
        next_sons=next_gen(:,j);
        next_sons(find(next_sons==0))=[];
        if ~isempty(next_sons)
            for k=(1:length(next_sons))
                index_s=find(order==next_sons(k));
                line([X_all(index_p),X_all(index_s)],[Y_all(index_p),Y_all(index_s)])
            end
        end
    end
end
end