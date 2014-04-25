function [ ] = GenWturbPlots(ConfigFile)
% GenWturbPlots generate plots showing both w var and w skew with dual x-axes

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Pdir = Config.PlotDir;

    Fsize = 25;

    % PlotDefs
    %
    % Format --> {
    %   input file prefix
    %   input file var name for w variance
    %   input file var name for w skew
    %   config plotset name
    %   x axis (variance) spec
    %   x axis (variance) spec
    %   y axis spec
    %   output file name
    %
    %   Axis spec:
    %     { [min max] scale label }
    PlotDefs = {
      {
      'turb_stats'
      'w-w_st_circ_all'
      'w-w-w_st_circ_all'
      'CO_UP_C0400_G10M5'
      '(a) C400'
      { [-5 5] 1e3  '<w\primew\prime> (m^2 s^-^2 x 10^-^3)' }
      { [-5 5] 1e3  '<w\primew\primew\prime> (m^3 s^-^3 x 10^-^3)' }
      { [ 0 4] 1e-3 'Height (km)' }
      'WturbStCirc_CO_C0400_G10M5_TALL.jpg'
      }

      };

    for iplot = 1:length(PlotDefs)
      InFprefix   = PlotDefs{iplot}{1};
      InWvarName  = PlotDefs{iplot}{2};
      InWskewName = PlotDefs{iplot}{3};
      PlotSetName = PlotDefs{iplot}{4};
      Ptitle      = PlotDefs{iplot}{5};
      XvarLim     = PlotDefs{iplot}{6}{1};
      XvarScale   = PlotDefs{iplot}{6}{2};
      XvarLabel   = PlotDefs{iplot}{6}{3};
      XskewLim    = PlotDefs{iplot}{7}{1};
      XskewScale  = PlotDefs{iplot}{7}{2};
      XskewLabel  = PlotDefs{iplot}{7}{3};
      Ylim        = PlotDefs{iplot}{8}{1};
      Yscale      = PlotDefs{iplot}{8}{2};
      Ylabel      = PlotDefs{iplot}{8}{3};
      OutFname    = PlotDefs{iplot}{9};

      fprintf('**************************************************************************\n');
      fprintf('Generating W tubulence statistics plot:\n');
      fprintf('  Plot Set: %s\n', PlotSetName);

      % find the index corresponding to the plot set name
      ips = FindPlotSet(Config.PlotSets, PlotSetName);
      if (ips <= 0)
        fprintf('ERROR: Cannot find plot set "%s" in configuration file "%s": skipping this plot.', PlotSetName, ConfigFile);
        fprintf('\n');
        break;
      else
        fprintf('  Plot Set Index: %d\n', ips);
        fprintf('\n');
      end

      % Read in the profiles for w variance and w skew. Place each profile in a column of
      % the corresponding array for the subsequent plot command.
      Nc = Config.PlotSets(ips).Ncases;
      for icase = 1:Nc
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText{icase} = Config.PlotSets(ips).Cases(icase).Legend;
        LineColors{icase} = Config.PlotSets(ips).Cases(icase).Lcolor;
        LineStyles{icase} = Config.PlotSets(ips).Cases(icase).Lstyle;
        LineGscales(icase) = Config.PlotSets(ips).Cases(icase).Lgscale;

        InFile = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);
        fprintf('  Reading: %s (%s, %s)\n', InFile, InWvarName, InWskewName);

        if (icase == 1)
          Z = squeeze(hdf5read(InFile, 'z_coords'));
          Nz = length(Z);

          WVAR  = zeros([ Nz Nc ]);
          WSKEW = zeros([ Nz Nc ]);
        end
        WVAR(:,icase)  = squeeze(hdf5read(InFile, InWvarName));
        WSKEW(:,icase) = squeeze(hdf5read(InFile, InWskewName));
      end
      fprintf('\n');

      WVAR = WVAR .* XvarScale;
      WSKEW = WSKEW .* XskewScale;
      Z = Z .* Yscale;

      % Create the plot
      %   Draw WVAR with solid lines and WSKEW with dashed lines
      %   Use two x axes - WVAR on bottom, WSKEW on top
      %     this forces two y axes - WVAR on left, WSKEW on right

      Fig = figure;

      % First plot the variance. Get the font sizes, labels, title, etc. all placed before
      % copying the axis to build the skew plot. This will get the axes positioned correctly
      % before copying them for the skew plot.
      plot(WVAR, Z, 'LineStyle', '-', 'LineWidth', 2);
      xlim(XvarLim);
      ylim(Ylim);

      set(gca, 'LineWidth', 2);
      set(gca, 'TickLength', [ 0.025 0.025 ]);

      set(gca, 'FontSize', Fsize);

      % create title
      T = title(Ptitle);
      set(T, 'Units', 'Normalized');
      set(T, 'HorizontalAlignment', 'Left');

      xlabel(XvarLabel);
      ylabel(Ylabel);

      % Shrink plot area so that double x-axis will fit
      Ppos = get(gca, 'Position'); % position of plot area
      Ppos(1) = Ppos(1) * 1.00;
      Ppos(2) = Ppos(2) * 1.00;
      Ppos(3) = Ppos(3) * 0.90;
      Ppos(4) = Ppos(4) * 0.80;
      set(gca, 'Position', Ppos);

      % move title over to the left side
      Tpos = get(T, 'Position');
      Tpos(1) = 0; % line up with left edge of plot area
      Tpos(2) = 1.25;  % move up above axis to make room for skew axis
      set(T, 'Position', Tpos);


      % Now put in skew axes and lines
      AxisPos = get(gca, 'Position');
      % Color = none makes the plot area transparent so var plots can be seen
      SkewAxis = axes('Position', AxisPos, 'XaxisLocation', 'top', ...
         'YaxisLocation', 'right', 'Color', 'none');
      line(WSKEW, Z, 'Parent', SkewAxis, 'LineStyle', '--', 'LineWidth', 2);

      set(SkewAxis, 'FontSize', Fsize);
      set(SkewAxis, 'Xlim', XskewLim, 'Ylim', Ylim);
      set(get(SkewAxis,'XLabel'), 'String', XskewLabel, 'FontSize', Fsize);
      set(get(SkewAxis,'YLabel'), 'String', Ylabel, 'FontSize', Fsize);
      
      OutFile = sprintf('%s/%s', Pdir, OutFname);
      fprintf('  Writing: %s\n', OutFile);
      saveas(Fig, OutFile);
      close(Fig); 
      fprintf('\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FindPlotSet()
%
% This function will return the index into the given list of plot sets corresponding
% to the given plot set name. If the plot set cannot be found then a zero is returned.
%
function [ Ind ] = FindPlotSet(Psets, PsName)

  Ind = 0;
  for i = 1:length(Psets)
    if (strcmp(Psets(i).Name, PsName))
      Ind = i;
      break;
    end
  end
end


%%%     % These lists are organized as 2D cell arrays. Each row is one complete spec for one variable.
%%%     % Syntax for rows:
%%%     %    { 'input file prefix' 'input file var name' 'input var  number' 'output var name' }
%%%     VarList = {
%%%       % means
%%%       { 'turb_mmts_w'               'turb_mmts_w'       1 'w'                }
%%%       { 'turb_cov_w_theta'          'turb_cov_theta'    2 'theta'            }
%%%       { 'turb_mmts_theta_e'         'turb_mmts_theta_e' 1 'theta_e'          }
%%%       { 'turb_cov_w_theta_v'        'turb_cov_theta_v'  2 'theta_v'          }
%%%       { 'turb_cov_w_vapor'          'turb_cov_vapor'    2 'vapor'            }
%%%       { 'turb_cov_w_speed'          'turb_cov_speed'    2 'speed'            }
%%%     
%%%       { 'turb_mmts_w_ud0p10'        'turb_mmts_w'       1 'w_ud0p10'         }
%%%       { 'turb_cov_w_theta_ud0p10'   'turb_cov_theta'    2 'theta_ud0p10'     }
%%%       { 'turb_mmts_theta_e_ud0p10'  'turb_mmts_theta_e' 1 'theta_e_ud0p10'   }
%%%       { 'turb_cov_w_theta_v_ud0p10' 'turb_cov_theta_v'  2 'theta_v_ud0p10'   }
%%%       { 'turb_cov_w_vapor_ud0p10'   'turb_cov_vapor'    2 'vapor_ud0p10'     }
%%%       { 'turb_cov_w_speed_ud0p10'   'turb_cov_speed'    2 'speed_ud0p10'     }
%%% 
%%%       { 'turb_mmts_w_st_circ'       'turb_mmts_w'       1 'w_st_circ'        }
%%%       { 'turb_mmts_theta_e_st_circ' 'turb_mmts_theta_e' 1 'theta_e_st_circ'  }
%%% 
%%%       { 'turb_mmts_w_cv_circ'       'turb_mmts_w'       1 'w_cv_circ'        }
%%%       { 'turb_mmts_theta_e_cv_circ' 'turb_mmts_theta_e' 1 'theta_e_cv_circ'  }
%%% 
%%%       % fluxes (covariances)
%%%       { 'turb_cov_w_theta'          'turb_cov_theta'    3 'w-theta'          }
%%%       { 'turb_cov_w_theta_v'        'turb_cov_theta_v'  3 'w-theta_v'        }
%%%       { 'turb_cov_w_vapor'          'turb_cov_vapor'    3 'w-vapor'          }
%%%       { 'turb_cov_w_speed'          'turb_cov_speed'    3 'w-speed'          }
%%%       
%%%       { 'turb_cov_w_theta_ud0p10'   'turb_cov_theta'    3 'w-theta_ud0p10'   }
%%%       { 'turb_cov_w_theta_v_ud0p10' 'turb_cov_theta_v'  3 'w-theta_v_ud0p10' }
%%%       { 'turb_cov_w_vapor_ud0p10'   'turb_cov_vapor'    3 'w-vapor_ud0p10'   }
%%%       { 'turb_cov_w_speed_ud0p10'   'turb_cov_speed'    3 'w-speed_ud0p10'   }
%%% 
%%%       % variances
%%%       { 'turb_mmts_w'               'turb_mmts_w'       2 'w-w'              }
%%% 
%%%       { 'turb_mmts_w_ud0p10'        'turb_mmts_w'       2 'w-w_ud0p10'       }
%%% 
%%%       { 'turb_mmts_w_st_circ'       'turb_mmts_w'       2 'w-w_st_circ'      }
%%% 
%%%       { 'turb_mmts_w_cv_circ'       'turb_mmts_w'       2 'w-w_cv_circ'      }
%%% 
%%%       % skews
%%%       { 'turb_mmts_w'               'turb_mmts_w'       3 'w-w-w'            }
%%% 
%%%       { 'turb_mmts_w_ud0p10'        'turb_mmts_w'       3 'w-w-w_ud0p10'     }
%%% 
%%%       { 'turb_mmts_w_st_circ'       'turb_mmts_w'       3 'w-w-w_st_circ'    }
%%% 
%%%       { 'turb_mmts_w_cv_circ'       'turb_mmts_w'       3 'w-w-w_cv_circ'    }
%%%       
%%%       };
%%% 
%%%     % do one hour time average at beginning, middle and end of sampling period
%%%     % also do time average over entire sampling inerval
%%%     TSstart = 12;
%%%     TSend = 13;
%%%     TMstart = 23.5;
%%%     TMend = 24.5;
%%%     TEstart = 35;
%%%     TEend = 36;
%%% 
%%%     TAstart = 12;
%%%     TAend = 36;
%%% 
%%%     TotalN = 158404; % excluding borders --> 398 * 398
%%% 
%%%     for icase = 1:length(Config.Cases)
%%%       Case = Config.Cases(icase).Cname;
%%%       OutFname = sprintf('%s/turb_stats_%s.h5', Ddir, Case);
%%% 
%%%       fprintf('***************************************************************\n');
%%%       fprintf('Generating moment/flux profiles:\n');
%%%       fprintf('  Case: %s\n', Case);
%%%       fprintf('  Output file: %s\n', OutFname);
%%%       fprintf('\n');
%%% 
%%%       % The vars and coords are written into the file in append mode so write 
%%%       % the file here without append mode in order to create a new file each time
%%%       % this script is run.
%%%       hdf5write(OutFname, 'header', Case);
%%% 
%%%       for ivar = 1:length(VarList)
%%%         InFprefix = VarList{ivar}{1};
%%%         InVname   = VarList{ivar}{2};
%%%         InOrder   = VarList{ivar}{3};
%%%         OutVname  = VarList{ivar}{4};
%%% 
%%%         InFname = sprintf('%s/%s_%s.h5', Tdir, InFprefix, Case);
%%% 
%%%         fprintf('  Input file: %s\n', InFname);
%%%         fprintf('    Var name: %s\n', InVname);
%%%         fprintf('    Var order: %d\n', InOrder);
%%%         fprintf('\n');
%%%         
%%%         % the number of points is always the 4th index of the first
%%%         % dimension
%%%         TURB_DATA = squeeze(hdf5read(InFname, InVname));
%%%         SUM = squeeze(TURB_DATA(InOrder,:,:,:));
%%%         N   = squeeze(TURB_DATA(4,:,:,:));
%%%         Z   = squeeze(hdf5read(InFname, 'z_coords'));
%%%         T   = squeeze(hdf5read(InFname, 't_coords')) / 3600;   % hr
%%%         
%%%         Nz = length(Z);
%%%         Nt = length(T);
%%% 
%%%         % find the indices of the start, mid, end and "all" time intervals
%%%         TS1 = find(T >= TSstart, 1, 'first');
%%%         TS2 = find(T <= TSend, 1, 'last');
%%% 
%%%         TM1 = find(T >= TMstart, 1, 'first');
%%%         TM2 = find(T <= TMend, 1, 'last');
%%% 
%%%         TE1 = find(T >= TEstart, 1, 'first');
%%%         TE2 = find(T <= TEend, 1, 'last');
%%% 
%%%         TA1 = find(T >= TAstart, 1, 'first');
%%%         TA2 = find(T <= TAend, 1, 'last');
%%% 
%%%         % N and SUM are organized as (z,t)
%%%         % Create the entire time series profiles, then average over time to
%%%         % get the start, middle, end and all profiles.
%%%         PROF = SUM ./ N;
%%%         
%%%         PROF_S = squeeze(nanmean(PROF(:,TS1:TS2),2));
%%%         PROF_M = squeeze(nanmean(PROF(:,TM1:TM2),2));
%%%         PROF_E = squeeze(nanmean(PROF(:,TE1:TE2),2));
%%%         PROF_A = squeeze(nanmean(PROF(:,TA1:TA2),2));
%%% 
%%%         % Generate a fraction statistic - the ratio of number of points selected 
%%%         % to total number of points in domain
%%%         PROF_FRAC = N ./ TotalN;
%%%         PROF_FRAC_S = squeeze(nanmean(PROF_FRAC(:,TS1:TS2), 2));
%%%         PROF_FRAC_M = squeeze(nanmean(PROF_FRAC(:,TM1:TM2), 2));
%%%         PROF_FRAC_E = squeeze(nanmean(PROF_FRAC(:,TE1:TE2), 2));
%%%         PROF_FRAC_A = squeeze(nanmean(PROF_FRAC(:,TA1:TA2), 2));
%%% 
%%%         % Write out data - put in dummy x, y and t coordinates
%%%         Xdummy = 1;
%%%         Ydummy = 1;
%%%         Tdummy = 1;
%%% 
%%%         OutVar = reshape(PROF_S, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_start', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_M, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_mid', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_E, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_end', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_A, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_all', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_FRAC_S, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_frac_start', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_FRAC_M, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_frac_mid', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_FRAC_E, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_frac_end', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         OutVar = reshape(PROF_FRAC_A, [ 1 1 Nz 1 ]);
%%%         OutVarName = sprintf('%s_frac_all', OutVname);
%%%         fprintf('  Writing var: %s\n', OutVarName);
%%%         hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');
%%% 
%%%         fprintf('\n');
%%%       end
%%% 
%%%       % all vars are written out to the file, now write out the coordinates
%%%       hdf5write(OutFname, 'x_coords', Xdummy, 'WriteMode', 'append');
%%%       hdf5write(OutFname, 'y_coords', Ydummy, 'WriteMode', 'append');
%%%       hdf5write(OutFname, 'z_coords', Z,      'WriteMode', 'append');
%%%       hdf5write(OutFname, 't_coords', Tdummy, 'WriteMode', 'append');
