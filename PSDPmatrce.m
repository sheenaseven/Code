function  Z=PSDPmatrce(positive_data,negtive_data)
%positive_data,negative_data are char array
[m_pos,n_pos]=size(positive_data);
pos_score=zeros(441,n_pos-1);  %441*(n_pos-1)æÿ’Û£¨”√”⁄¥Ê∑≈PSDPæÿ’Û
V=[];
k=1;
for j=1:n_pos-1
    Seq=positive_data(:,j:j+1);
    F=dipetide_composition(Seq);
    F=sum(F);
    pos_score(:,j)=F';
    Seq=[];
end
 [m_neg,n_neg]=size(negtive_data);
 for i=1:10
    p=randperm(m_neg)';
    negtive_data=negtive_data(p,:);
    neg=negtive_data(1:m_pos,:);
    [m,n]=size(neg);
    for j=1:n-1
        Seq=neg(:,j:j+1);
        F=dipetide_composition(Seq);
        F=sum(F);
        neg_score(i).mat(:,j)=F';
        Seq=[];
    end
 end
 for j=1:n_pos-1
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