function [ ] = EqMeasHdf5()

  SimList = {
%    'RCE_3km_1mom'
%    'RCE_3km_1mom_db'
    'RCE_3km_1mom_db_udef'
    'RCE_3km_1mom_db_rlongup'
%    'RCE_3km_1mom_dm'
%    'RCE_3km_2mom_db'
    'RCE_3km_2mom_db_udef'
    'RCE_3km_2mom_db_rlongup'
%    'RCE_3km_2mom_dm'
%    'RCE_3km_2mom_dm_lrz'
    };
  Nsims = length(SimList);
  
  for i = 1:Nsims
    Sim = SimList{i};
  
    % Read in matlab file
    InFile = sprintf('%s/Lfiles/EnergyBalance.mat', Sim);
    fprintf('Loading: %s\n', InFile);
  
    % After loading file:
    %   meanIntRadCool is the net radiative flux time series
    %   meanSFlux is the surface flux time series
    load(InFile);
  
    OutFile = sprintf('eq_meas_%s.h5', Sim);
    RfluxVname = '/net_rad_flux';
    SfluxVname = '/sfc_flux';
    IradFluxVname = '/int_rad_flux';
  
    ShfVname = '/shf';
    LhfVname = '/lhf';
  
    LwupTopVname = '/lwup_top';
    SwupTopVname = '/swup_top';
    SwdnTopVname = '/swdn_top';
  
    LwupSfcVname = '/lwup_sfc';
    LwdnSfcVname = '/lwdn_sfc';
    SwupSfcVname = '/swup_sfc';
    SwdnSfcVname = '/swdn_sfc';
  
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end
  
    WriteDataset(OutFile, RfluxVname, meanNetRad);
    WriteDataset(OutFile, SfluxVname, meanSFlux);
    WriteDataset(OutFile, IradFluxVname, meanIntRadCool);
  
    WriteDataset(OutFile, ShfVname, meanSHF);
    WriteDataset(OutFile, LhfVname, meanLHF);
  
    WriteDataset(OutFile, LwupTopVname, meanLWupTop);
    WriteDataset(OutFile, SwupTopVname, meanSWupTop);
    WriteDataset(OutFile, SwdnTopVname, meanSWdnTop);
  
    WriteDataset(OutFile, LwupSfcVname, meanLWupSfc);
    WriteDataset(OutFile, LwdnSfcVname, meanLWdnSfc);
    WriteDataset(OutFile, SwupSfcVname, meanSWupSfc);
    WriteDataset(OutFile, SwdnSfcVname, meanSWdnSfc);
  
    fprintf('\n');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = WriteDataset(OutFname, OutVname, OutVar)

  fprintf('Writing: %s (%s)\n', OutFname, OutVname);
  h5create(OutFname, OutVname, length(OutVar));
  h5write(OutFname, OutVname, OutVar);

end
