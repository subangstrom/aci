function [xydist,h]=position_analysis(fitresult,AImage,range,edge)
psize=length(fitresult);
[ImageX, ImageY]=size(AImage);

for i=1:1:psize
    for j=1:1:psize
        if(length(fitresult{j})<7 || length(fitresult{i})<7)
            continue;
        end
        if(fitresult{i}(6)<1+edge || fitresult{i}(6)>ImageX-edge || fitresult{i}(5)<1+edge ||fitresult{i}(5)>ImageY-edge)
            continue;
        end
        if(fitresult{j}(6)<1+edge || fitresult{j}(6)>ImageX-edge || fitresult{j}(5)<1+edge ||fitresult{j}(5)>ImageY-edge)
            continue;
        end
        xydist(i,j)=sqrt((fitresult{i}(5)-fitresult{j}(5))^2+(fitresult{i}(6)-fitresult{j}(6))^2);
    end
end
[sizex,sizey]=size(xydist);
xydist=reshape(xydist,1,sizex*sizey);
r=0:0.1:range;
h=hist(xydist,r);
h(1)=0;
h(length(h))=0;
end
