function  X=PSAAPcode(Z,Seq)

model=['A'    'R'    'N'    'D'    'C'    'Q'    'E'  ...
    'G'    'H'    'I'    'L'    'K'    'M'    'F' ...  
    'P'    'S'    'T'    'W'    'Y'    'V'];

% [a,b]=size(Seq);
%  Seq=Seq(:,[1:fix(b/2),fix(b/2)+2:end]);

[m_seq,n_seq]=size(Seq);
X=[];
for i=1:m_seq
    for j=1:n_seq
        z=find(Seq(i,j)==model);
        if(length(z)==0)
          disp( Seq(i,j));
        else
            X(i,j)=Z(z,j);
        end
    end
end
        











