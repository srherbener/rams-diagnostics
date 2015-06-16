function [ ] = GenBarGraphFiles()
% GenBarGraphFiles generate data for bar graphs

    Ddir = 'DIAGS';

    VarSets = {
      {
        'Average cloud optical thickness'
        {
          { 'DIAGS/avg_cot_<CASE>.h5' '/cot_tavg'         '/cot_avg_<CSET>'         }
          { 'DIAGS/avg_cot_<CASE>.h5' '/cot_all_cld_tavg' '/cot_all_cld_avg_<CSET>' }
        }
        'DIAGS/bgraph_cot.h5'
      }

      };
    Nvsets = length(VarSets);

    CaseSets = {
      {
        's293'
        {
          'z.atex.ccn0050.sst293'
          'z.atex.ccn0100.sst293'
          'z.atex.ccn0200.sst293'
          'z.atex.ccn0400.sst293'
          'z.atex.ccn0800.sst293'
          'z.atex.ccn1600.sst293'
        }
      }

      {
        's298'
        {
          'z.atex.ccn0050.sst298'
          'z.atex.ccn0100.sst298'
          'z.atex.ccn0200.sst298'
          'z.atex.ccn0400.sst298'
          'z.atex.ccn0800.sst298'
          'z.atex.ccn1600.sst298'
        }
      }

      };
    Ncsets = length(CaseSets);

    for ivset = Nvsets
      Descrip = VarSets{ivset}{1};
      VarList = VarSets{ivset}{2};
      OutFile = VarSets{ivset}{3};

      Nvars = length(VarList);

      fprintf('*************************************************************************\n');
      fprintf('Generating bar graph data: %s\n', Descrip);
      fprintf('\n');

      % Remove the output file if it exists so that the writes below will not collide
      % with datasets left in an old version of the output file.
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      for ivar = 1:Nvars
        InFtemplate  = VarList{ivar}{1};
        InVname      = VarList{ivar}{2};
        OutVtemplate = VarList{ivar}{3};

        for icset = 1:Ncsets
          CsetName = CaseSets{icset}{1};
          CaseList = CaseSets{icset}{2};
   
          Ncases = length(CaseList);

          Vdata = zeros([ 1 Ncases ]);
          for icase = 1:Ncases
            Case         = CaseList{icase};

            InFile = regexprep(InFtemplate, '<CASE>', Case);
            fprintf('  Reading: %s (%s)\n', InFile, InVname);
            Vdata(icase) = squeeze(h5read(InFile, InVname));
          end
          fprintf('\n');

          OutVname = regexprep(OutVtemplate, '<CSET>', CsetName);
          fprintf('  Writing: %s (%s)\n', OutFile, OutVname);

          h5create(OutFile, OutVname, size(Vdata));
          h5write (OutFile, OutVname, Vdata);

          fprintf('\n');
        end
      end
    end

      % dump out the averages plus the coordinate values into the output file

%      hdf5write(OutFile, 'Averages', OutAvgs);
%      hdf5write(OutFile, 'Npoints',  OutNpts,  'WriteMode', 'append');
%
%      hdf5write(OutFile, 'VarNames', VarList,   'WriteMode', 'append');
%      hdf5write(OutFile, 'SST',      SstVals,  'WriteMode', 'append');
%      hdf5write(OutFile, 'CCN',      CcnVals,  'WriteMode', 'append');
%      hdf5write(OutFile, 'GCCN',     GccnVals, 'WriteMode', 'append');

end
