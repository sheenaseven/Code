function M=KSpace(seqs)
%对序列矩阵按space间隔的氨基酸对进行编码
%seqs 序列矩阵，行数为序列条数,短肽数据
%A=[1,2];A=[1]
%space 氨基酸对间隔数
amino =['A'    'C'    'D'    'E'    'F'    'G'    'H'  'I'    'K'  'L'    'M'    'N'    'P'    'Q'   'R'    'S'    'T'    'V'    'W'    'Y'   'X' ];

%amino =['R'    'G'    'I'    'F'    'P'    'S'    'T'  'Y'    'V' ]; 'L'    'M'    'N'    'P'    'Q'   'R'    'S'    'T'    'V'    'W'    'Y'   'X'
% matrix_code=zeros(length(amino),length(amino),size(seqs,1));
M=zeros(size(seqs,1),length(amino)*length(amino),size(seqs,2)-1);  %size(seqs,1)返回矩阵的行
for space=0:size(seqs,2)-2;%所有space的情况
 matrix_code=zeros(length(amino),length(amino),size(seqs,1));
   for j = 1:size(seqs,1)
     seq_singal = seqs(j,:);  %取seqs的每一行短肽ACKEF或者EFCKA
     for i=1:size(seqs,2)-space-1  %对matrix_code进行基于space的两个氨基酸编码
        a1=find(amino==seq_singal(i));
        a2=find(amino==seq_singal(i+space+1));
        matrix_code(a1,a2,j)=matrix_code(a1,a2,j)+1/(length(amino)*length(amino));
        
     end    
   end
 for m=1:size(seqs,1)  
    sub_code(m,:) = reshape(matrix_code(:,:,m)',1,length(amino)*length(amino));  %按照AA，AC，AE等等对应成一行441列
 end
    M(:,:,space+1)=sub_code(:,:);%将sub_code按短肽个数排成几行
end

