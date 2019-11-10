%% Parameters

lambda = 1
resolution = 0.01
Tilt_x = -10*lambda
Tilt_y = -10*lambda
colors = gray;
%% Polar Space

xh = -1:resolution:1;
[x,y] = meshgrid(xh,xh);
r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);
in = double(gt(1,r));
in(in==0) = NaN;

%% Zernike Polynomials

Z_00 = 1;
Z_11 = 2.*r.*cos(theta);
Z_1n1 = 2.*r.*sin(theta);
Z_20 = sqrt(3).*(2.*(r.^2) - 1);
Z_22 = sqrt(6).*(r.^2).*cos(2.*theta);
Z_2n2 = sqrt(6).*(r.^2).*sin(2.*theta);
Z_31 = 2*sqrt(2).*(3.*(r.^3) - 2.*r).*cos(theta);
Z_3n1 = 2*sqrt(2).*(3.*(r.^3) - 2.*r).*sin(theta);
Z_33 = 2*sqrt(2).*(r.^3).*cos(3.*theta);
Z_3n3 = 2*sqrt(2).*(r.^3).*sin(3.*theta);
Z_40 = sqrt(5).*(6.*(r.^4) - 6.*(r.^2)+1);


%% Wavefronts

%Reference wave
E_ref = exp((2*pi*1i/lambda).*(Tilt_x.*x + Tilt_y.*y));
% Test wave
E_r = in.*exp((2*pi*1i/lambda));
I = (abs(E_r+E_ref).^2)./max((abs(E_r+E_ref).^2),[],'all');
% Interference Pattern of spherical wave
image(xh,xh,I,'CDataMapping','scaled')
set(gca,'visible','off')
colormap(gray)
shading interp 


%% Aberration

abb_options = {Z_00, Z_11, Z_1n1, Z_20, Z_22, Z_2n2, Z_31, Z_3n1, Z_33, Z_3n3, Z_40};
% If you want to decide on specific coefficients
abb_c = lambda*[0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

% if you need random coefficients
%{
abb_c = rand(1,11)
%}
abb = 0;
for i = 1:11
    abb = abb + abb_c(i)*cell2mat(abb_options(i));
end


E_abb = in.*exp((2*pi*1i/lambda).*abb);
surfl(x,y,in.*abb)
set(gca,'visible','off')
colormap(gray)
shading interp 
imshow(in.*abb,[])
I_abb = (abs(E_abb+E_ref).^2)./max((abs(E_abb+E_ref).^2),[],'all');
%image(xh,xh,I_abb,'CDataMapping','scaled')

%colormap(colors)
%shading interp 
%colorbar

%% RMS

in = double(gt(1,r));
in(in==0) = NaN;
abb = in.*abb;
M = nanmean(abb,'all')
RMS = (sqrt(nanmean(abb.^2,'all') - nanmean(abb,'all').^2))/lambda
RMS_predicted = sqrt(sum(abb_c(2:11).^2))/lambda

image(xh,xh,I_abb,'CDataMapping','scaled')
text(-0.7,1.1,sprintf('Predicted RMS = %f', RMS_predicted),'FontSize',14)
text(-0.7,1.2,sprintf('Calculated RMS = %f', RMS),'FontSize',14)
colormap(colors)
shading interp 
set(gca,'visible','off')