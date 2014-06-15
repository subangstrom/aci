
function [fitresult,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center,C,D,StartPoint,h_area]=find_atomic_columns(raw_image,sigma,threshold,max_peak_num,sign,style,area_threshold,initial_values,fit_shift,verbose)

%sign=1 for adf
%sign=-1 for abf
%style=1 to use the correlated image for peak finding
%style=2 to use the original image for peak finding

warning('off','optimlib:lsqncommon:SwitchToLargeScale');

[ImageX, ImageY]=size(raw_image);
%find the best template used for searching peaks (atomic columns)
if(sigma==0)
    template_search_range=10;
    maxC=zeros(template_search_range,1);
    for i=1:1:template_search_range
        template=sign*fspecial('gaussian',20,i);
        C=normxcorr2(template,raw_image);
        maxC(i,1)=max(max(C));
    end
    sigma=find(maxC==max(maxC));
end
fprintf('gaussian template sigma=%d\n',sigma);

% use the default threshold to detect peaks
if(threshold==0)
    threshold=0.75;
end

for i=1:1:1
    template=sign*fspecial('gaussian',19,sigma);
    C=normxcorr2(template,raw_image);%ser.data{1,i});
    maxC=max(max(C));
    D=C(10:10+ImageX-1,10:10+ImageY-1);  % this might be a different way to get the peak position
    %D=ImageSum;
    C(C<threshold*maxC)=0;
    C(C>threshold*maxC)=1;
    image_peak=raw_image.*C(10:10+ImageX-1,10:10+ImageY-1);
    if verbose ==1
        figure;
        colormap(gray);
        imagesc(image_peak);
        axis image;
    end
    CC=bwconncomp(image_peak);
    areas_in_pixels = cellfun(@length, CC.PixelIdxList);
    fprintf('total peaks found: %d areas range from %d to %d\n',CC.NumObjects,min(areas_in_pixels),max(areas_in_pixels));
    fprintf('peaks with areas larger than 10 pixels: %d\n',length(areas_in_pixels(areas_in_pixels>10)));
    fprintf('peaks with areas larger than 20 pixels: %d\n',length(areas_in_pixels(areas_in_pixels>20)));
    
    area_size=[1:1:1000];
    h_area=hist(areas_in_pixels,area_size);
    
    if verbose ==1
        figure;
        plot(h_area);
    end
    
    % use the centroid first, see if we can improve the result later using
    % fitting
    centroid = regionprops(CC, 'centroid');
    
    %area_threshold=250;
    expand=0; %expand the area to include more pixels
    
    %     xpos=zeros(CC.NumObjects,1);
    %     ypos=zeros(CC.NumObjects,1);
    %     intensity=zeros(CC.NumObjects,1);
    image1=zeros(ImageX,ImageY);
    image2=zeros(ImageX,ImageY);
    indexj=1;
    for j=1:1:min(CC.NumObjects,max_peak_num)
        if(areas_in_pixels(j)<area_threshold(1) || areas_in_pixels(j)>area_threshold(2) || areas_in_pixels(j)>1000 || areas_in_pixels(j)<10)
            continue;
        end
        peak_index=CC.PixelIdxList{j};
        peak_x=mod(peak_index,ImageX);
        peak_x(peak_x==0)=ImageX;
        peak_y=floor((peak_index-1)/ImageX)+1;
        
%         peak_x_min=max(min(peak_x)-expand,1);
%         peak_x_max=min(max(peak_x)+expand,ImageX);
%         peak_y_min=max(min(peak_y)-expand,1);
%         peak_y_max=min(max(peak_y)+expand,ImageY);
%         fy=peak_x_min:1:peak_x_max;
%         fx=peak_y_min:1:peak_y_max;
        %peak_data=ser.data{1,1}(peak_x_min:peak_x_max,peak_y_min:peak_y_max);
        %[peak_x_min,peak_x_max,peak_y_min,peak_y_max]
        fx=peak_x;
        fy=peak_y;
        fz=peak_x;
        for pixels=1:1:length(peak_index)
            if style==1
                fz(pixels)=sign*D(fx(pixels),fy(pixels));
            end
            if style==2
                fz(pixels)=sign*raw_image(fx(pixels),fy(pixels));
            end
        end
%         if style==1
%             peak_data=sign*D(peak_x_min:peak_x_max,peak_y_min:peak_y_max);
%         end
%         if style==2
%             peak_data=sign*raw_image(peak_x_min:peak_x_max,peak_y_min:peak_y_max);
%         end

        
%         fz=peak_data;
        if ~iscell(initial_values)
            [fitresult{indexj},zfit{indexj}, fiterr(indexj,1:7), zerr{indexj}, resnorm(indexj), rr(indexj),StartPoint{indexj}] = fmgaussfit_improved(fy,fx,fz,0,fit_shift);
        else
            [fitresult{indexj},zfit{indexj}, fiterr(indexj,1:7), zerr{indexj}, resnorm(indexj), rr(indexj),StartPoint{indexj}] = fmgaussfit_improved(fy,fx,fz,initial_values{indexj},fit_shift);
        end
        mass_center(indexj,1:2)=centroid(j).Centroid;
        object_index(indexj)=j;
        indexj=indexj+1;
        
        for pixels=1:1:length(peak_index)
            image1(fx(pixels),fy(pixels))=fz(pixels);
            image2(fx(pixels),fy(pixels))=zfit{indexj-1}(pixels);
        end
        oimage{indexj}=fz;
        fprintf('find peak number %d in frame %d\n',j,i);
    end
%     figure;
%     imagesc(image1);
%     figure;
%     imagesc(image2);
end
end