function [ ] = PlotInitialAerosols(ConfigFile)
% PlotInitialAerosols Plot the aerosol profiles used in the ATEX simulations

Config = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

% Using the vertical profile mode of aerosol initialization (ICLOUD = 6)
CCN = [ 50 100 200 400 800 1200 1600 ];
GCCN = [ 1e-5 1e-4 1e-2 1 ];

ClegText = {
    '50/cc'
    '100/cc'
    '200/cc'
    '400/cc'
    '800/cc'
    '1200/cc'
    '1600/cc'
    };

GlegText = {
    '10^-^5/cc'
    '10^-^4/cc'
    '10^-^2/cc'
    '1/cc'
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

plot(C_PROFS, Z, 'LineWidth', 3);
set(gca, 'FontSize', Fsize);
xlabel('N_C_C_N (#/cc)');
ylabel('Height (m)');

legend(ClegText, 'Location', 'NorthEast', 'FontSize', 20);
legend boxoff;

OutFile = sprintf('%s/InitialCcnProfile.jpg', Pdir);
fprintf('Writing file: %s\n', OutFile);
saveas(Fig, OutFile);

close(Fig);

% GCCN profiles
Fig = figure;

semilogx(G_PROFS, Z, 'LineWidth', 3);
set(gca, 'FontSize', Fsize);
set(gca, 'Xtick', [ 1e-5 1e-3 1e-1 ]);
xlabel('N_G_C_C_N (#/cc)');
ylabel('Height (m)');
xlim([ 1e-6 1 ]);
ylim([ 0 7000 ]);

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
