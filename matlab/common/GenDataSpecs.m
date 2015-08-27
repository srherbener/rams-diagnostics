function [ DataSpecs, LegText, DSokay ] = GenDataSpecs(Config, FigCase, i_panel, i_pset, Indent)
% GenDataSpecs function to generate data specfifications for plots

  UndefVal = Config.UndefVal;
  Smooth   = Config.FigPanels(i_panel).Smooth;
  Flength  = Config.FigPanels(i_panel).Flength;
  Ndsets   = Config.PlotSets(i_pset).Ndsets;

  % Figure out which data sets (X, Y, Z) are to be smoothed
  % Smooth is a string that specifies which data sets to operate on.
  %   'xyz' -> smooth x y and z sets
  %   'xz'  -> smooth x and z, but not y
  %   'y'   -> smooth only y
  Xsmooth = ~isempty(regexp(Smooth,'x'));
  Ysmooth = ~isempty(regexp(Smooth,'y'));
  Zsmooth = ~isempty(regexp(Smooth,'z'));

  DSokay = 1;

  for i_dset = 1:Ndsets
    i_pdata = Config.PlotSets(i_pset).DataSets(i_dset).PDnum;
    if (i_pdata == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotData on PlotSet number %d, Line number %d\n', i_panel, i_pset, i_dset)
      DSokay = 0;
      continue;
    end
 
    % Line properties
    DataSpecs(i_dset).Lwidth  = Config.PlotSets(i_pset).DataSets(i_dset).Lwidth;
    DataSpecs(i_dset).Lcolor  = Config.PlotSets(i_pset).DataSets(i_dset).Lcolor;
    DataSpecs(i_dset).Lstyle  = Config.PlotSets(i_pset).DataSets(i_dset).Lstyle;
    DataSpecs(i_dset).Lgscale = Config.PlotSets(i_pset).DataSets(i_dset).Lgscale;

    % Legend text for this line
    LegText{i_dset} = Config.PlotSets(i_pset).DataSets(i_dset).Legend;

    % Line specific case
    % The line specific case takes precedence over the figure case
    LineCase = Config.PlotSets(i_pset).DataSets(i_dset).Case;
    if (strcmp(LineCase, 'none'))
      Case = FigCase;
    else
      Case = LineCase;
    end

    % Nvars indicates how many variables are specified for this
    % PlotData set. Assume the ordering:
    %   Vars(1) -> X data
    %   Vars(2) -> Y data
    %   Vars(3) -> Z data
    Nvars = Config.PlotData(i_pdata).Nvars;

    if (Nvars >= 1)
      DataSpecs(i_dset).Xdata = ReadProcessVarData(Config.PlotData(i_pdata).Vars(1), Case, Xsmooth, Flength, UndefVal, Indent);
    end
    if (Nvars >= 2)
      DataSpecs(i_dset).Ydata = ReadProcessVarData(Config.PlotData(i_pdata).Vars(2), Case, Ysmooth, Flength, UndefVal, Indent);
    end
    if (Nvars >= 3)
      DataSpecs(i_dset).Zdata = ReadProcessVarData(Config.PlotData(i_pdata).Vars(3), Case, Zsmooth, Flength, UndefVal, Indent);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ReadProcessVarData()
%
% This function will pull the specs for reading data for a variable
% out of the config structure, read in the variable, apply the
% appropriate processing and output the result.
%
function [ VarData ] = ReadProcessVarData( Var, Case, Smooth, Flength, UndefVal, Indent )

    File   = regexprep(Var.File, '<CASE>', Case);
    
    % Read in var data
    if (strcmp(Var.Vname, 'dummy'))
      % use the pattern in File to generate VarData
      fprintf('%sUsing pattern: %s\n', Indent, File);
      VarData = eval(File);
    else
      % read in VarData from the file
      fprintf('%sReading: %s (%s)\n', Indent, File, Var.Vname);
      VarData = ReadSelectHdf5(File, Var.Vname, Var.Select);
    end

    % Apply specifice processing
    %    linear scaling
    %    undefined values set to nan
    VarData = (VarData .* Var.Scale) + Var.Offset;
    VarData(VarData == UndefVal) = nan;
    if (Smooth)
      VarData = SmoothFillTseries(VarData, length(VarData), Flength);
    end

end
