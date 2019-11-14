%% Properties
resolution = 0.005;
wavelength = 1;
amplitude = 1;
xh = -1:resolution:1;

%% Define Waves
g = geometry(resolution);
testph = height_error(g,{'R'},g.circa);
%testph.plot_height(g);
t = wave(wavelength,testph.distribution,amplitude,g); %test wave definition
refph = height_error(g,{'T',10,10},1);
r = wave(wavelength,refph.distribution,amplitude,g); %reference wave definition

%% Superposition
sp_real = PSI(t,r,0,g,'centre');
sp_real.c_nm
sp_pure = PSI(t,r,3,g,'centre');
sp_pure.c_nm
sp_pure.plot_interference(g);
%sp.plot_phase(g)






