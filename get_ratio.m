function [ratio,ratio_matrix,col_int_matrix]=get_ratio(fitresult,mini,col_int,a1,a2)
psize=length(fitresult);
[sx,sy]=size(mini);
ratio=zeros(psize,1);
ratio_matrix=zeros(sx,sy);
col_int_matrix=zeros(sx,sy);

for i=1+a1:1:sx-a1
    for j=1+a2:1:sy-a2
        if(mini(i,j)==0)
            continue;
        end
        p1=mini(i,j);
        
        p2=mini(i+a1,j+a2);
        p3=mini(i+a1,j-a2);
        p4=mini(i-a1,j+a2);
        p5=mini(i-a1,j-a2);
        
        if(p2==0 || p3==0 || p4==0 || p5==0)
            continue;
        end
        
        ratio(p1)=col_int(p1)*4/(col_int(p2)+col_int(p3)+col_int(p4)+col_int(p5));
        ratio_matrix(i,j)=ratio(p1);
        col_int_matrix(i,j)=col_int(p1);
    end
end
end