function  Z=PSAAPmatrce(positive_data,negtive_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%统计正的样本中每个位置各种氨基酸的频率
[m_pos,n_pos]=size(positive_data);
pos_score=zeros(21,n_pos);  %与PSDP的主要区别处pos_score=zeros(441,n_pos-1);
V=[];
k=1;
for i=1:n_pos
    pos_score(:,i)=amino_acid_composition(positive_data(:,i)')';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m_neg,n_neg]=size(negtive_data);
for i=1:10
    p=randperm(m_neg)';
    negtive_data=negtive_data(p,:);
    neg=negtive_data(1:m_pos,:);
    [m,n]=size(neg);
    for j=1:n
        neg_score(i).mat(:,j)=amino_acid_composition(neg(:,j)')';
    end
end

    for j=1:n_pos
       for i=1:10
        V=[V  neg_score(i).mat(:,j)];
       end  
        neg_mean(:,k)=mean(V')';
        neg_std(:,k)=std(V')';
        V=[];
        k=k+1;
   end
[m_score,n_score]=size(pos_score);
for i=1:m_score
    for j=1:n_score
        if neg_std(i,j)==0
            Z(i,j)=0;
        else
            Z(i,j)=(pos_score(i,j)-neg_mean(i,j))/neg_std(i,j);
        end
    end
end
%%%%%%%%%%%%%%%%%中间字符一样,所以去除的结果
% b=size(Z,2);
% Z(:,(fix(b/2)+1))=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%