function code_seq=code_amino(seqs,kernel,times,space)
%function code_seq=code_amino(seqs,positive_data,negtive_data,times,space)
%给多条序列编码
%seqs 序列矩阵，行数为序列条数
%positive_data 正类点序列矩阵，行数为序列条数
%negtive_data  负类点序列矩阵，行数为序列条数
%times 负类点序列矩阵随机抽取的次数,一般取10
%space 序列间隔数
for j=1:size(seqs,1)
    seq=seqs(j,:);
    code_index_seq = [];
    for i=space
        code_index_seq=[code_index_seq;code_amino_pair(seq,i)];
    end
    code_index_seq=code_index_seq';
    code_seq(j,:)=sum(code_index_seq.*kernel);
end
end


