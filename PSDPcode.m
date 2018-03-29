function  X=PSDPcode(Z,Seq)
[m_seq,n_seq]=size(Seq);
seq=[];
for i=1:m_seq
   
        for l=1:n_seq-1
        seq(l,:)=Seq(i,l:l+1);
        end
       F=dipetide_composition(seq);
       seq=[];
       [m,n]=size(F);
      
       for j=1:m
          z=find(F(j,:)==1);
           X(i,j)=Z(z,j);
       end
    
end
       