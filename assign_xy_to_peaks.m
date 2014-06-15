function [d,proj_acc,projx,peak_index,row_map,col_map,mini,row_stat,col_stat,coord_angle]=assign_xy_to_peaks(AImage,fitresult,cutoff1,A1,A2,count_threshold,tolerance)
psize=length(fitresult);
[ImageX, ImageY]=size(AImage);
if(A1==A2)
    d=[1:0.2:360];
else
    d=[A1-5:0.1:A1+5,A2-5:0.1:A2+5];
end
%d=[1:1:180];
dsize=length(d);
diag=floor(sqrt(ImageX*ImageX+ImageY*ImageY))+1;
proj=zeros(dsize,diag);
proj_acc=zeros(dsize,diag);
edge=30;

for j=1:1:dsize
    s=sind(d(j));
    c=cosd(d(j));
    spoint=floor(min([0,ImageX*cosd(d(j)),ImageX*cosd(d(j))+ImageY*sind(d(j)),ImageY*sind(d(j))]));
    epoint=floor(max([0,ImageX*cosd(d(j)),ImageX*cosd(d(j))+ImageY*sind(d(j)),ImageY*sind(d(j))]));
    for i=1:1:psize
        if(length(fitresult{i})<7)
            continue;
        end
        if(fitresult{i}(6)<1+edge || fitresult{i}(6)>ImageX-edge || fitresult{i}(5)<1+edge ||fitresult{i}(5)>ImageY-edge)
            continue;
        end
        x=fitresult{i}(6);
        y=fitresult{i}(5);
        a1=abs(s*x+c*y);
        a3=sqrt(x*x+y*y);
        
        a2=floor(x*cosd(d(j))+y*sind(d(j)));
        proj(j,a2-spoint+1)=1;
        proj_acc(j,a2-spoint+1)=proj_acc(j,a2-spoint+1)+1;
    end
end

%nor_I=ones(ImageX,ImageY);
%Iaccu=radon(nor_I,d);

%proj_acc=proj_acc./


%projx=sum(proj');
projx=std(proj_acc');
if cutoff1==0
    cutoff1=300;
end

projx1=projx(1:dsize/2);
projx2=projx(dsize/2+1:dsize);

cutoff1=mean(projx1)-std(projx1);
projmask1=blur(projx1,1,21);

projmask1(projmask1<cutoff1)=-100*cutoff1;
[ymax1,imax1,ymin1,imin1] = extrema(projmask1);

cutoff2=mean(projx2)-std(projx2);
projmask2=blur(projx2,1,21);

projmask2(projmask2<cutoff2)=-100*cutoff2;
[ymax2,imax2,ymin2,imin2] = extrema(projmask2);

imin(1)=imax1(1);
imin(2)=imax2(1)+dsize/2;
% if(imin(2)<imin(1))
%     temp=imin(1);
%     imin(1)=imin(2);
%     imin(2)=temp;
% end
for i=1:1:length(imin)
    fprintf('find image aligned along %f degree at index %d\n',d(imin(i)),imin(i));
end

coord_angle(1)=d(imin(1));
coord_angle(2)=d(imin(2));
axis1=proj_acc(imin(1),:);
axis2=proj_acc(imin(2),:);

%axis1(axis1<1)=0;
%axis2(axis2<1)=0;

axis1=conv(axis1, [1 2 1], 'same');
axis2=conv(axis2, [1 2 1], 'same');

%axis1(axis1<count_threshold)=0;
%axis2(axis2<count_threshold)=0;





% figure,plot(axis1);
% figure,plot(axis2);

axis1_C=bwlabeln(axis1);
axis2_C=bwlabeln(axis2);

axisj(1,:)=axis1_C;
axisj(2,:)=axis2_C;

peak_index=zeros(psize,2);
row_map=zeros(ImageX, ImageY);
col_map=zeros(ImageX, ImageY);

rc=100;
row_stat=zeros(rc,rc);
col_stat=zeros(rc,rc);

for j=1:1:2
    s=sind(d(imin(j)));
    c=cosd(d(imin(j)));
    spoint=floor(min([0,ImageX*c,ImageX*c+ImageY*s,ImageY*s]));
    epoint=floor(max([0,ImageX*c,ImageX*c+ImageY*s,ImageY*s]));
    for i=1:1:psize
        if(length(fitresult{i})<7)
            continue;
        end
        if(fitresult{i}(6)<1+edge || fitresult{i}(6)>ImageX-edge || fitresult{i}(5)<1+edge ||fitresult{i}(5)>ImageY-edge)
            continue;
        end
        x=fitresult{i}(6);
        y=fitresult{i}(5);
        a1=abs(s*x+c*y);
        a3=sqrt(x*x+y*y);
        
        a2=floor(x*c+y*s);
        px=floor(fitresult{i}(6));
        py=floor(fitresult{i}(5));
        
        peak_index(i,j)=axisj(j,a2-spoint+1);
        if j==1
            row_map(px,py)=max(axisj(j,a2-spoint+1-tolerance:a2-spoint+1+tolerance));
        end
        if j==2
            col_map(px,py)=max(axisj(j,a2-spoint+1-tolerance:a2-spoint+1+tolerance));
        end
    end
end

for i=1:1:psize
    if(length(fitresult{i})<7)
        continue;
    end
    if(fitresult{i}(6)<1+edge || fitresult{i}(6)>ImageX-edge || fitresult{i}(5)<1+edge ||fitresult{i}(5)>ImageY-edge)
        continue;
    end
    if(peak_index(i,1)==0 || peak_index(i,2)==0)
        continue;
    end
    mini(peak_index(i,1),peak_index(i,2))=i;
end

[minix,miniy]=size(mini);
row_stat=zeros(minix,miniy);
col_stat=zeros(minix,miniy);

%x value for each peak on the x-axis
for i=1:1:minix
    for j=1:1:miniy
        index=mini(i,j);
        if(index==0)
            continue;
        end
        s=sind(d(imin(1)));
        c=cosd(d(imin(1)));
        x=fitresult{index}(6);
        y=fitresult{index}(5);
        a1=abs(s*x+c*y);
        a3=sqrt(x*x+y*y);
        a2=sqrt(a3*a3-a1*a1);
        row_stat(i,j)=a2;
    end
end

% y value for each peak on the y-axis
for i=1:1:minix
    for j=1:1:miniy
        index=mini(i,j);
        if(index==0)
            continue;
        end
        s=sind(d(imin(2)));
        c=cosd(d(imin(2)));
        x=fitresult{index}(6);
        y=fitresult{index}(5);
        a1=abs(s*x+c*y);
        a3=sqrt(x*x+y*y);
        a2=sqrt(a3*a3-a1*a1);
        col_stat(i,j)=a2;
    end
end

end
