%% Properties
resolution = 0.001;
wavelength = 1;
amplitude = 1;
xh = -1:resolution:1;

%% Define Waves
g = geometry(resolution);
testph = height_error(g,{'Z','Z_31'},g.circa);
%testph.plot_height(g);
t = wave(wavelength,testph.distribution,amplitude,g); %test wave definition
refph = height_error(g,{'T',10,10},1);
r = wave(wavelength,refph.distribution,amplitude,g); %reference wave definition

%% Superposition
sp = PSI(t,r,'example',g,'centre');
sp.plot_interference(g);
sp.plot_phase(g)






