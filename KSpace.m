function M=KSpace(seqs)
%�����о���space����İ�����Խ��б���
%seqs ���о�������Ϊ��������,��������
%A=[1,2];A=[1]
%space ������Լ����
amino =['A'    'C'    'D'    'E'    'F'    'G'    'H'  'I'    'K'  'L'    'M'    'N'    'P'    'Q'   'R'    'S'    'T'    'V'    'W'    'Y'   'X' ];

%amino =['R'    'G'    'I'    'F'    'P'    'S'    'T'  'Y'    'V' ]; 'L'    'M'    'N'    'P'    'Q'   'R'    'S'    'T'    'V'    'W'    'Y'   'X'
% matrix_code=zeros(length(amino),length(amino),size(seqs,1));
M=zeros(size(seqs,1),length(amino)*length(amino),size(seqs,2)-1);  %size(seqs,1)���ؾ������
for space=0:size(seqs,2)-2;%����space�����
 matrix_code=zeros(length(amino),length(amino),size(seqs,1));
   for j = 1:size(seqs,1)
     seq_singal = seqs(j,:);  %ȡseqs��ÿһ�ж���ACKEF����EFCKA
     for i=1:size(seqs,2)-space-1  %��matrix_code���л���space���������������
        a1=find(amino==seq_singal(i));
        a2=find(amino==seq_singal(i+space+1));
        matrix_code(a1,a2,j)=matrix_code(a1,a2,j)+1/(length(amino)*length(amino));
        
     end    
   end
 for m=1:size(seqs,1)  
    sub_code(m,:) = reshape(matrix_code(:,:,m)',1,length(amino)*length(amino));  %����AA��AC��AE�ȵȶ�Ӧ��һ��441��
 end
    M(:,:,space+1)=sub_code(:,:);%��sub_code�����ĸ����ųɼ���
end

