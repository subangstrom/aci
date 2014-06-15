function [Iproj,Iavg,Istd]=project_image_RD(I,threshold,d)

tic

[ImageX, ImageY]=size(I);
%d=[-5:0.1:5,85:0.1:95];
%d=[-90+angle:1:90+angle];
dsize=length(d);
nor_I=ones(ImageX,ImageY);
Iproj=radon(I,d);
Iaccu=radon(nor_I,d);
Iproj=Iproj./Iaccu;

for j=1:1:dsize
    
    temp=Iproj(:,j);
    temp_c=Iaccu(:,j);
    temp(temp_c<threshold)=0;
    temp=temp(temp~=0);
    
    Iavg(j)=mean(temp);
    Istd(j)=std(temp);
end

Iprojx=sum(Iproj');
toc
Istd=blur(Istd);

end

function b=blur(a)
b=a;
for i=1+1:1:length(a)-1
    b(i)=a(i-1)*0.1+a(i)+a(i+1)*0.1;
end
end