function code_seq=code_amino(seqs,kernel,times,space)
%function code_seq=code_amino(seqs,positive_data,negtive_data,times,space)
%���������б���
%seqs ���о�������Ϊ��������
%positive_data ��������о�������Ϊ��������
%negtive_data  ��������о�������Ϊ��������
%times ��������о��������ȡ�Ĵ���,һ��ȡ10
%space ���м����
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


