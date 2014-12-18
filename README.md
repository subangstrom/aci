Atom Column Indexing (ACI)
===

This folder contains matlab files that load an example RevSTEM image and then fit atom columns in STEM image using 2D Gaussian distribution and then index the atom columns to a matrix representation.

The matlab file "ACI_flow.m" contains the sections that can run in sequence by using the command "Run section" in Matlab software. The same result should be expected as shown in "html/ACI_flow.html".

<html><head>
      <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <title></title><meta name="generator" content="MATLAB 8.1"><link rel="schema.DC" href="http://purl.org/dc/elements/1.1/"><meta name="DC.date" content="2014-12-11"><meta name="DC.source" content="ACI_flow.m"><style type="text/css">
html,body,div,span,applet,object,iframe,h1,h2,h3,h4,h5,h6,p,blockquote,pre,a,abbr,acronym,address,big,cite,code,del,dfn,em,font,img,ins,kbd,q,s,samp,small,strike,strong,sub,sup,tt,var,b,u,i,center,dl,dt,dd,ol,ul,li,fieldset,form,label,legend,table,caption,tbody,tfoot,thead,tr,th,td{margin:0;padding:0;border:0;outline:0;font-size:100%;vertical-align:baseline;background:transparent}body{line-height:1}ol,ul{list-style:none}blockquote,q{quotes:none}blockquote:before,blockquote:after,q:before,q:after{content:'';content:none}:focus{outine:0}ins{text-decoration:none}del{text-decoration:line-through}table{border-collapse:collapse;border-spacing:0}

html { min-height:100%; margin-bottom:1px; }
html body { height:100%; margin:0px; font-family:Arial, Helvetica, sans-serif; font-size:10px; color:#000; line-height:140%; background:#fff none; overflow-y:scroll; }
html body td { vertical-align:top; text-align:left; }

h1 { padding:0px; margin:0px 0px 25px; font-family:Arial, Helvetica, sans-serif; font-size:1.5em; color:#d55000; line-height:100%; font-weight:normal; }
h2 { padding:0px; margin:0px 0px 8px; font-family:Arial, Helvetica, sans-serif; font-size:1.2em; color:#000; font-weight:bold; line-height:140%; border-bottom:1px solid #d6d4d4; display:block; }
h3 { padding:0px; margin:0px 0px 5px; font-family:Arial, Helvetica, sans-serif; font-size:1.1em; color:#000; font-weight:bold; line-height:140%; }

a { color:#005fce; text-decoration:none; }
a:hover { color:#005fce; text-decoration:underline; }
a:visited { color:#004aa0; text-decoration:none; }

p { padding:0px; margin:0px 0px 20px; }
img { padding:0px; margin:0px 0px 20px; border:none; }
p img, pre img, tt img, li img { margin-bottom:0px; } 

ul { padding:0px; margin:0px 0px 20px 23px; list-style:square; }
ul li { padding:0px; margin:0px 0px 7px 0px; }
ul li ul { padding:5px 0px 0px; margin:0px 0px 7px 23px; }
ul li ol li { list-style:decimal; }
ol { padding:0px; margin:0px 0px 20px 0px; list-style:decimal; }
ol li { padding:0px; margin:0px 0px 7px 23px; list-style-type:decimal; }
ol li ol { padding:5px 0px 0px; margin:0px 0px 7px 0px; }
ol li ol li { list-style-type:lower-alpha; }
ol li ul { padding-top:7px; }
ol li ul li { list-style:square; }

.content { font-size:1.2em; line-height:140%; padding: 20px; }

pre, tt, code { font-size:12px; }
pre { margin:0px 0px 20px; }
pre.error { color:red; }
pre.codeinput { padding:10px; border:1px solid #d3d3d3; background:#f7f7f7; }
pre.codeoutput { padding:10px 11px; margin:0px 0px 20px; color:#4c4c4c; }

@media print { pre.codeinput, pre.codeoutput { word-wrap:break-word; width:100%; } }

span.keyword { color:#0000FF }
span.comment { color:#228B22 }
span.string { color:#A020F0 }
span.untermstring { color:#B20000 }
span.syscmd { color:#B28C00 }

.footer { width:auto; padding:10px 0px; margin:25px 0px 0px; border-top:1px dotted #878787; font-size:0.8em; line-height:140%; font-style:italic; color:#878787; text-align:left; float:none; }
.footer p { margin:0px; }
.footer a { color:#878787; }
.footer a:hover { color:#878787; text-decoration:underline; }
.footer a:visited { color:#878787; }

table th { padding:7px 5px; text-align:left; vertical-align:middle; border: 1px solid #d6d4d4; font-weight:bold; }
table td { padding:7px 5px; text-align:left; vertical-align:top; border:1px solid #d6d4d4; }





  </style></head><body><div class="content"><h2>Contents</h2><div><ul><li><a href="#1">load the RevSTEM image. The gray-scale image is loaded to the variable 'ImageSum_R'</a></li><li><a href="#2">display the RevSTEM image</a></li><li><a href="#3">find the atom column locations using normalized cross correlation data</a></li><li><a href="#4">find the atom column locations using experimental RevSTEM data</a></li><li><a href="#5">distance histogram calculated from the fitting result</a></li><li><a href="#6">calculate PSD for the example image</a></li><li><a href="#7">we can then use the point-PSD to find exactly locations of peaks at roughly 85 degrees and -5 degrees (175 degrees in the plot)</a></li><li><a href="#8">note the output from matlab says 'find image aligned a long xxx degree at index yyy', use the index yyy to plot the projected profile</a></li><li><a href="#9">the matrix representation is stored in 'mini_E'</a></li></ul></div><h2>load the RevSTEM image. The gray-scale image is loaded to the variable 'ImageSum_R'<a name="1"></a></h2><p>the RevSTEM image has not been filtered note that our code can also work on regular STEM images the 'serReader.m' script can be used to read a converntiaonl STEM image acquired by TIA</p><pre class="codeinput">load <span class="string">example.mat</span>
</pre><h2>display the RevSTEM image<a name="2"></a></h2><pre class="codeinput">figure;
imagesc(ImageSum_R);
axis <span class="string">image</span>;
colormap(gray);
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_01.png" alt=""> <h2>find the atom column locations using normalized cross correlation data<a name="3"></a></h2><p>function [fitresult,oimage,zfit, fiterr, zerr, resnorm,... rr,image1,image2,object_index,mass_center,C,D,StartPoint,h_area]=... find_atomic_columns(raw_image,sigma,threshold,max_peak_num,sign,style,... area_threshold,initial_values,fit_shift,verbose);</p><pre class="codeinput"><span class="comment">% the main output 'fitresult_N' has a cell structure,</span>
<span class="comment">% each cell is a 1x7 array [amp, ang, sx, sy, xo, yo, zo] containing the peak fitting result</span>
<span class="comment">% amp: amplitude</span>
<span class="comment">% ang: rotation angle of the two main axes</span>
<span class="comment">% sx: sigma along the first main axis</span>
<span class="comment">% sy: sigma along the second main axis</span>
<span class="comment">% xo: x coordinate of the peak center</span>
<span class="comment">% yo: y coordinate of the peak center</span>
<span class="comment">% zo: background intensity</span>

<span class="comment">% Meaning of inputs of 'find_atom_clolumns'</span>
<span class="comment">% raw_image: The input STEM image</span>
<span class="comment">% sigma: sigma of the gaussian distribution for normalized cross-correlation,</span>
<span class="comment">% when the number is 0, the program decides the best sigma for ncc</span>
<span class="comment">% threshold: threshold to separate the atom columns</span>
<span class="comment">% max_peak_num: the limit of peak numbers</span>
<span class="comment">% sign: 1: find the peaks; 2: find the valleys (for ABF)</span>
<span class="comment">% style: 1: use ncc data to fit the peak; 2: use experimental data</span>
<span class="comment">% area_threshold: only areas within the area_threshold range are used for fitting</span>
<span class="comment">% initial_values: starting values for peak fitting, can be set as 0</span>
<span class="comment">% fit_shift: for future use, 0</span>
<span class="comment">% verbose: show the ncc threshold map and area size histogram</span>
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center_N]=<span class="keyword">...</span>
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,1,[100 400],0,0,1);
</pre><pre class="codeoutput">gaussian template sigma=8
total peaks found: 974 areas range from 154 to 37563
peaks with areas larger than 10 pixels: 974
peaks with areas larger than 20 pixels: 974
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_02.png" alt=""> <img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_03.png" alt=""> <h2>find the atom column locations using experimental RevSTEM data<a name="4"></a></h2><p>style is set to 2 to use experimental data for fitting starting points are set to be fitresult_N</p><pre class="codeinput">[fitresult_E,oimage_E,zfit_E, fiterr_E, zerr_E, resnorm_E, rr_E,image1_E,image2_E,startpoint_E,mass_center_E]=<span class="keyword">...</span>
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,2,[100 400],fitresult_N,0,0);
</pre><pre class="codeoutput">gaussian template sigma=8
total peaks found: 974 areas range from 154 to 37563
peaks with areas larger than 10 pixels: 974
peaks with areas larger than 20 pixels: 974
</pre><h2>distance histogram calculated from the fitting result<a name="5"></a></h2><pre class="codeinput">figure;
[xydist_E,h_E]=position_analysis(fitresult_E,ImageSum_R,300,30);
plot(h_E);
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_04.png" alt=""> <h2>calculate PSD for the example image<a name="6"></a></h2><pre class="codeinput">figure;
[Iproj,Iavg,Istd]=project_image_RD(ImageSum_R,100,-90:1:90);
plot(Istd);
</pre><pre class="codeoutput">Elapsed time is 5.352220 seconds.
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_05.png" alt=""> <h2>we can then use the point-PSD to find exactly locations of peaks at roughly 85 degrees and -5 degrees (175 degrees in the plot)<a name="7"></a></h2><pre class="codeinput">figure;
[d_E,proj_acc_E,projx_E,peak_index_E,row_map_E,col_map_E,mini_E,row_stat_E,col_stat_E,coord_angle_E]<span class="keyword">...</span>
    =assign_xy_to_peaks(ImageSum_R,fitresult_E,200,-5,85,1,1);
plot(projx_E);
</pre><pre class="codeoutput">find image aligned along -6.100000 degree at index 40
find image aligned along 84.100000 degree at index 143
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_06.png" alt=""> <h2>note the output from matlab says 'find image aligned a long xxx degree at index yyy', use the index yyy to plot the projected profile<a name="8"></a></h2><pre class="codeinput">figure;
subplot(2,1,1); plot(proj_acc_E(40,:));
subplot(2,1,2); plot(proj_acc_E(143,:));
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_07.png" alt=""> <h2>the matrix representation is stored in 'mini_E'<a name="9"></a></h2><p>mini_E is a 2D matrix with each node containing the index of the peak fitting result in fitresult_R</p><pre class="codeinput">figure;
imagesc(mini_E&gt;0);
daspect([1 sqrt(2) 1]);
colormap(gray);
</pre><img vspace="5" hspace="5" src="https://raw.githubusercontent.com/subangstrom/aci/master/html/ACI_flow_08.png" alt=""> <p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R2013a</a><br></p></div><!--
##### SOURCE BEGIN #####
%% load the RevSTEM image. The gray-scale image is loaded to the variable 'ImageSum_R'
% the RevSTEM image has not been filtered
% note that our code can also work on regular STEM images
% the 'serReader.m' script can be used to read a converntiaonl STEM image
% acquired by TIA

load example.mat

%% display the RevSTEM image

figure;
imagesc(ImageSum_R);
axis image;
colormap(gray);

%% find the atom column locations using normalized cross correlation data
% function [fitresult,oimage,zfit, fiterr, zerr, resnorm,...
% rr,image1,image2,object_index,mass_center,C,D,StartPoint,h_area]=...
% find_atomic_columns(raw_image,sigma,threshold,max_peak_num,sign,style,...
% area_threshold,initial_values,fit_shift,verbose);

% the main output 'fitresult_N' has a cell structure,
% each cell is a 1x7 array [amp, ang, sx, sy, xo, yo, zo] containing the peak fitting result
% amp: amplitude
% ang: rotation angle of the two main axes
% sx: sigma along the first main axis
% sy: sigma along the second main axis
% xo: x coordinate of the peak center
% yo: y coordinate of the peak center
% zo: background intensity

% Meaning of inputs of 'find_atom_clolumns'
% raw_image: The input STEM image
% sigma: sigma of the gaussian distribution for normalized cross-correlation,
% when the number is 0, the program decides the best sigma for ncc
% threshold: threshold to separate the atom columns
% max_peak_num: the limit of peak numbers
% sign: 1: find the peaks; 2: find the valleys (for ABF)
% style: 1: use ncc data to fit the peak; 2: use experimental data
% area_threshold: only areas within the area_threshold range are used for fitting 
% initial_values: starting values for peak fitting, can be set as 0
% fit_shift: for future use, 0
% verbose: show the ncc threshold map and area size histogram
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center_N]=...
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,1,[100 400],0,0,1);

%% find the atom column locations using experimental RevSTEM data
% style is set to 2 to use experimental data for fitting
% starting points are set to be fitresult_N
[fitresult_E,oimage_E,zfit_E, fiterr_E, zerr_E, resnorm_E, rr_E,image1_E,image2_E,startpoint_E,mass_center_E]=...
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,2,[100 400],fitresult_N,0,0);

%% distance histogram calculated from the fitting result
figure;
[xydist_E,h_E]=position_analysis(fitresult_E,ImageSum_R,300,30);
plot(h_E);

%% calculate PSD for the example image
figure;
[Iproj,Iavg,Istd]=project_image_RD(ImageSum_R,100,-90:1:90);
plot(Istd);

%% we can then use the point-PSD to find exactly locations of peaks at roughly 85 degrees and -5 degrees (175 degrees in the plot) 
figure;
[d_E,proj_acc_E,projx_E,peak_index_E,row_map_E,col_map_E,mini_E,row_stat_E,col_stat_E,coord_angle_E]...
    =assign_xy_to_peaks(ImageSum_R,fitresult_E,200,-5,85,1,1);
plot(projx_E);

%% note the output from matlab says 'find image aligned a long xxx degree at index yyy', use the index yyy to plot the projected profile

figure;
subplot(2,1,1); plot(proj_acc_E(40,:));
subplot(2,1,2); plot(proj_acc_E(143,:));

%% the matrix representation is stored in 'mini_E'
% mini_E is a 2D matrix with each node containing the index of the peak
% fitting result in fitresult_R

figure;
imagesc(mini_E>0);
daspect([1 sqrt(2) 1]);
colormap(gray);

##### SOURCE END #####
--></body></html>
