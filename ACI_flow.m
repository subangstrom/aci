%% load the RevSTEM image. The data is loaded to the variable 'ImageSum_R'
load example.mat

%% display the RevSTEM image
figure;
imagesc(ImageSum_R);
axis image;
colormap(gray);

%% find the atom column locations using normalized cross correlation data
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center_N]=...
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,1,[100 400],0,0,0);

%% find the atom column locations using experimental RevSTEM data
[fitresult_E,oimage_E,zfit_E, fiterr_E, zerr_E, resnorm_E, rr_E,image1_E,image2_E,startpoint_E,mass_center_E]=...
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,2,[100 400],fitresult_N,0,0);

%% distance histogram calculated from the fitting
[xydist_E,h_E]=position_analysis(fitresult_E,ImageSum_R,300,30);
plot(h_E);

%% calculate PSD for the example image
[Iproj,Iavg,Istd]=project_image_RD(ImageSum_R,100,-90:1:90);
plot(Istd);

%% we can then use the point-PSD to find exactly locations of peaks at roughly 85 degrees and -5 degrees (175 degrees in the plot) 
[d_E,proj_acc_E,projx_E,peak_index_E,row_map_E,col_map_E,mini_E,row_stat_E,col_stat_E,coord_angle_E]...
    =assign_xy_to_peaks(ImageSum_R,fitresult_E,200,-5,85,1,1);
plot(projx_E);

%% note the output from matlab says 'find image aligned a long xxx degree at index yyy', use the index yyy to plot the projected profile
subplot(2,1,1); plot(proj_acc_E(40,:));
subplot(2,1,2); plot(proj_acc_E(143,:));

%% the matrix representation is stored in 'mini_E'
figure;
imagesc(mini_E>0);
daspect([1 sqrt(2) 1]);
colormap(gray);
