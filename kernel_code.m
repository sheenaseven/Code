function kernel=kernel_code(positive_data,negtive_data,times,space)
%建立核矩阵
%positive_data 正类点序列矩阵，行数为序列条数
%negtive_data  负类点序列矩阵，行数为序列条数
%times 负类点序列矩阵随机抽取的次数
%space 氨基酸对间隔数
code_pos=[];
for i=space
    code_pos=[code_pos;code_amino_pair(positive_data,i)];
end

for j=1:times
    rn=randperm(size(negtive_data,1));
    negtive_data_use=negtive_data(rn(1:size(positive_data,1)),:);
    code_neg_temp=[];
    for i=space       
        code_neg_temp=[code_neg_temp;code_amino_pair(negtive_data_use,i)];
    end
    code_neg(:,:,j)=code_neg_temp;
end

code_mean=zeros(size(code_pos));
for j=1:times
    code_mean=code_mean+code_neg(:,:,j);
end
code_mean=code_mean/times;

code_std=zeros(size(code_pos));
for j=1:times
    code_std=code_std+(code_neg(:,:,j)-code_mean).^2;
end
code=(code_pos-code_mean)./code_std; 
code(code_std==0)=0;
kernel=code';
end