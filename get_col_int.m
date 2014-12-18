function [ col_int ] = get_col_int( AImage,fitresult,range,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
psize=length(fitresult);

col_int=zeros(psize,1);

[ImageX, ImageY]=size(AImage);

%edge=2;
for i=1:1:psize
    if(length(fitresult{i})<7)
        continue;
    end
    if(fitresult{i}(6)<1+2 || fitresult{i}(6)>ImageX-2 || fitresult{i}(5)<1+2 ||fitresult{i}(5)>ImageY-2)
        continue;
    end
    if(method==1)
        col_int(i)=intensity(AImage,fitresult{i}(6),fitresult{i}(5),range);
    end
    if(method==0)
    %as the experimental data is directly used for fitting, we can use the
    %fittng result to get the peak intensity
        col_int(i)=fitresult{i}(1)+fitresult{i}(7);
    end
end

end

function inte=intensity(AImage,i,j,range)
ir=range;
inte=0;
[ImageX, ImageY]=size(AImage);
for k=-ir:1:ir
    for l=-ir:1:ir
        if (k*k+l*l>range*range)
            continue;
        end
        indexx=i+k;
        indexy=j+l;
        if(indexx>1 && indexx<ImageX && indexy>1 && indexy<ImageY)
            inte=inte+cubic_interpolate(AImage,indexx,indexy);
        end
    end
end
end

function inte=cubic_interpolate(AImage,i,j)
A=AImage(floor(i)-2:floor(i)+2,floor(j)-2:floor(j)+2);
x=floor(i)-2:floor(i)+2;
y=floor(j)-2:floor(j)+2;
for k= 1:1:5
    z(k)=csapi(y,A(k,:),j);
end
inte=csapi(x,z,i);
end
