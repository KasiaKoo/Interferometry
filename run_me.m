%% Properties
resolution = 0.005;
wavelength = 1;
phase = 1;
amplitude = 1;
xh = -1:resolution:1;

%% Define Object

g = geometry(resolution);
t = wave(wavelength,phase,amplitude,g); %test wave definition
t = t.aberration('Z_31',g); %add aberation to test wave
r = wave(wavelength,phase,amplitude,g); %reference wave definition
r = r.tilt(5,5,g); %add tilt to reference wave

%% Superposition
figure()
intens=abs(t.front+r.front).^2;
imshow(intens, [])
t.plot_height(2,g)
rms = sqrt(nanmean(g.circa.*t.height_error(:).^2)-nanmean(g.circa.*t.height_error(:)));
coma = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*cos(g.theta).*g.circa;
c = nansum(g.circa.*t.height_error(:).*coma(:))/nansum(g.circ(:));




