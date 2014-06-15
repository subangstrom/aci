function [fitresult, zfit, fiterr, zerr, resnorm, rr,StartPoint] = fmgaussfit_improved(xx,yy,zz,initial_values,fit_shift)
% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface 
%   height. XX and YY are vectors or matrices defining the x and y 
%   components of a surface.
%          
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, zfit, fiterr, zerr, resnorm, rr] =
%       fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.

%   Copyright 2013, Nathan Orloff.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% link to the mathworks download page http://www.mathworks.com/matlabcentral/fileexchange/41938-fit-2d-gaussian-with-optimization-toolbox


%% Condition the data
xData=xx;
yData=yy;
xyData = {xData,yData};
zData=zz;

%% Set up the startpoint 
[amp, ind] = max(zData); % amp is the amplitude. 
if(initial_values==0)
    xo =  mean(xData);%xData(ind); % guess that it is at the maximum
    yo =  mean(yData);%yData(ind); % guess that it is at the maximum
    ang = 0; % angle in degrees.
    sy = 5;
    sx = 5;
    
else
    xo=initial_values(5);
    yo=initial_values(6);
    sx=initial_values(3)+2;
    sy=initial_values(4)+2;
    ang=initial_values(2);
end
zo = median(zData(:))-std(zData(:));
xmax = max(xData)+2; 
ymax = max(yData)+2; 
zmax = amp*2; % amp is the amplitude. 
xmin = min(xData)-2; 
ymin = min(yData)-2; 
zmin = min(zData)/2; % amp is the amplitude.

amp=amp-zo;

%% Set up fittype and options. 
Lower = [0, eps, 0, 0, xmin, ymin, zmin]; 
Upper = [Inf, 180, Inf, Inf, xmax, ymax, zmax]; % angles greater than 90 are redundant 
StartPoint = [amp, ang, sx, sy, xo, yo, zo];%[amp, sx, sxy, sy, xo, yo, zo];

tols = 1e-16;
options = optimset('Algorithm','levenberg-marquardt',...
    'Display','off',...
    'MaxFunEvals',5e2,...
    'MaxIter',5e2,...
    'TolX',tols,...
    'TolFun',tols,...
    'TolCon',tols ,...
    'UseParallel','always');

%% perform the fitting
[fitresult,resnorm,residual] = ...
    lsqcurvefit(@gaussian2D,StartPoint,xyData,zData,Lower,Upper,options);

sizexy=length(xyData{1});

offsetX = ones(sizexy,1)*fit_shift;
offsetY = ones(sizexy,1)*fit_shift;

offset = {offsetX+xyData{1}, offsetY+xyData{2}};

[fiterr, zfit, zerr] = gaussian2Duncert(fitresult,residual,offset);
rr = rsquared(zData, zfit, zerr);
zfit = reshape(zfit,size(zz));
zerr = reshape(zerr,size(zz));

end

function rr = rsquared(z,zf,ze)
% reduced chi-squared
dz = z-zf;
rr = 1/(numel(z)-8)*sum(dz.^2./ze.^2); % minus 8 because there are 7 fit parameters +1 (DOF)
end

function z = gaussian2D(par,xy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
end

function [dpar,zf,dzf] = gaussian2Duncert(par,resid,xy)
% get the confidence intervals
J = guassian2DJacobian(par,xy);
parci = nlparci(par,resid,'Jacobian',J);
dpar = (diff(parci,[],2)./2)';
[zf,dzf] = nlpredci(@gaussian2D,xy,par,resid,'Jacobian',J);
end

function J = guassian2DJacobian(par,xy)
% compute the jacobian
x = xy{1}; y = xy{2};
J(:,1) = exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2);
J(:,2) = -par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(3).^2 - (2.*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2);
J(:,3) = (2.*par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2)./par(3)^3;
J(:,4) = (2.*par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2)./par(4)^3;
J(:,5) = par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*cosd(par(2)).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))))./par(3).^2 - (2.*sind(par(2)).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2);
J(:,6) = par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*cosd(par(2)).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2 + (2.*sind(par(2)).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))))./par(3).^2);
J(:,7) = ones(size(x));
end