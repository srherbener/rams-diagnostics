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

Z = 50:100:4950;

Cnmin = 50;
Gnmin = 1e-5;

Nc = length(CCN);
Ng = length(GCCN);
Nz = length(Z);

C_PROFS = zeros(Nz, Nc);
G_PROFS = zeros(Nz, Ng);

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

OutFile = sprintf('%/InitialCcnProfile.jpg', Pdir);
saveas(Fig, OutFile);

close(Fig);

% GCCN profiles
Fig = figure;

plot(G_PROFS, Z, 'LineWidth', 3);

OutFile = sprintf('%/InitialGccnProfile.jpg', Pdir);
saveas(Fig, OutFile);

close(Fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenRamsProfile() 
%
% This function will emulate the code in RAMS that creates
% the initial vertical profiles for the aerosols.

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
  if ((k > Kmin) & (Z(k) <= Zmax))
    Prof(k) = max(Nmin, Nmax * (1 - (Z(k)/Zmax)));
  end
  if (Z(k) > Zmax)
    Prof(k) = Nmin;
  end
end

end
