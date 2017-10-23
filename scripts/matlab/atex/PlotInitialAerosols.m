function [ ] = PlotInitialAerosols(ConfigFile)
% PlotInitialAerosols Plot the aerosol profiles used in the ATEX simulations

Config = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

% Using the vertical profile mode of aerosol initialization (ICLOUD = 6)
CCN = [ 50 100 200 400 800 1600 ];
GCCN = [ 1e-4 1e-2 1 ];

ClegText = {
    'C50'
    'C100'
    'C200'
    'C400'
    'C800'
    'C1600'
    };

GlegText = {
    'G10M4'
    'G10M2'
    'G10M0'
    };

CcolorList = {
    'black'
    'green'
    'blue'
    'cyan'
    'red'
    'magenta'
    };

GcolorList = {
    'cyan'
    'blue'
    'magenta'
    };

Z = 50:100:4950;

Cnmin = 50;
Gnmin = 1e-5;

Nc = length(CCN);
Ng = length(GCCN);
Nz = length(Z);

C_PROFS = zeros(Nz, Nc);
G_PROFS = zeros(Nz, Ng);

Fsize = 35;

for i = 1:Nc
  C_PROFS(:,i) = GenRamsProfile(Z, Cnmin, CCN(i), 1, 4000);
end

for i = 1:Ng
  G_PROFS(:,i) = GenRamsProfile(Z, Gnmin, GCCN(i), 3, 4000);
end

% Make the plots

% CCN profiles
Fig = figure;

Hlines = plot(C_PROFS, Z./1000, 'LineWidth', 3);
for i = 1:Nc
  set(Hlines(i), 'Color', str2rgb(CcolorList{i}));
end
set(gca, 'FontSize', Fsize);
xlabel('N_a (# cm^-^3)');
ylabel('Height (km)');
ylim([ 0 5 ]);

legend(ClegText, 'Location', 'NorthEast', 'FontSize', 20);
legend boxoff;

OutFile = sprintf('%s/InitialCcnProfile.jpg', Pdir);
fprintf('Writing file: %s\n', OutFile);
saveas(Fig, OutFile);

close(Fig);

% GCCN profiles
Fig = figure;

Hlines = semilogx(G_PROFS, Z./1000, 'LineWidth', 3);
for i = 1:Ng
  set(Hlines(i), 'Color', str2rgb(GcolorList{i}));
end
set(gca, 'FontSize', Fsize);
set(gca, 'Xtick', [ 1e-5 1e-3 1e-1 ]);
xlabel('N_a (# cm^-^3)');
ylabel('Height (km)');
xlim([ 1e-6 1 ]);
ylim([ 0 7 ]);

legend(GlegText, 'Location', 'NorthWest', 'Orientation', 'horizontal', 'FontSize', 20);
legend boxoff;

OutFile = sprintf('%s/InitialGccnProfile.jpg', Pdir);
fprintf('Writing file: %s\n', OutFile);
saveas(Fig, OutFile);

close(Fig);

end


function [ Prof ] = GenRamsProfile(Z, Nmin, Nmax, Kmin, Zmax)
% GenRamsProfile emulate the code in RAMS that generates the vertical profile for aerosol concentration
%    Z    - vector containing z levels
%    Nmin - minimum number concentration
%    Nmax - maximum number concentration
%    Kmin - lower bound (index of Z vector) of ramp in the profile
%    Zmax - upper bound (height, units of Z) of ramp in the profile

Nz = length(Z);

Prof = zeros([ Nz 1 ]);
for k = 1:Nz
  if (k <= Kmin)
    Prof(k) = Nmax;
  end
  if ((k > Kmin) && (Z(k) <= Zmax))
    Prof(k) = max(Nmin, Nmax * (1 - (Z(k)/Zmax)));
  end
  if (Z(k) > Zmax)
    Prof(k) = Nmin;
  end
end

end
