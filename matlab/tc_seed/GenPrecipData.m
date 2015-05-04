function [] = GenPrecipData()
% GenPrecipData fucntion to generate time averaged precip rate 

  Ddir = 'DIAGS';
  
  % Make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
      mkdir(Ddir);
  end
  
  Vname = 'pcprate';
  Hvar = sprintf('/%s', Vname);

  CaseList = {
    'TCS_SD_C0100'
    'TCS_SD_C0500'
    'TCS_SD_C1000'
    'TCS_SD_C2000'
    };
  
  % Steady state time period (hours)
  Tstart = 120;
  Tend   = 140;
  
  for icase = 1:length(CaseList);
    Case = CaseList{icase};
  
    Hfile = sprintf('AzAveragedData/pcprate_%s.h5', Case);
    OutFile = sprintf('%s/%s_%s_SS.h5', Ddir, Vname, Case);
  
    fprintf('***********************************************************************\n');
    fprintf('Generating radar reflectivity data:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Input file: %s\n', Hfile);
    fprintf('  Output file: %s\n', OutFile);
    fprintf('\n');
  
    PR = squeeze(h5read(Hfile, Hvar));
    R  = squeeze(h5read(Hfile, '/x_coords')) ./ 1000; % km
    T  = squeeze(h5read(Hfile, '/t_coords')) ./ 3600; % h

    % Find the indices for the time interval
    T1 = find(T >= Tstart, 1, 'first');
    T2 = find(T <= Tend  , 1, 'last');

    % Form the time average. PR is (r,t)
    AVG_PR = squeeze(nanmean(PR(:,T1:T2), 2))'; % make a column vector
    Nr = length(AVG_PR);
  
    % write output file with attached dims so that grads can read the file
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    % AVG_PR is now (r)
  
    h5create(OutFile, Hvar, length(AVG_PR));
    h5write(OutFile, Hvar, AVG_PR);
  
    % write coordinate vars
    h5create(OutFile, '/radius', length(R));
    h5write(OutFile, '/radius', R);
  
    fprintf('\n');
  end
end
