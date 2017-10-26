function [] = FixSstFiles()

  % This routine will reshape the 2D data SEATF inside the SST files. This is done
  % since the data was created before the dimension swapping issue was repaired
  % in RAMS.

  SstPath = '/Volumes/tasman/sherbener/NOBAK/ACTIVE/simdata/TsDebby/SPIN_UP_SIM/MODEL';

  FileList = {
    'ssth-W-2006-08-09-000000-g1.h5'
    'ssth-W-2006-08-09-000000-g2.h5'
    'ssth-W-2006-08-09-000000-g3.h5'
    'ssth-W-2006-08-17-000000-g1.h5'
    'ssth-W-2006-08-17-000000-g2.h5'
    'ssth-W-2006-08-17-000000-g3.h5'
    'ssth-W-2006-08-25-000000-g1.h5'
    'ssth-W-2006-08-25-000000-g2.h5'
    'ssth-W-2006-08-25-000000-g3.h5'
    'ssth-W-2006-09-02-000000-g1.h5'
    'ssth-W-2006-09-02-000000-g2.h5'
    'ssth-W-2006-09-02-000000-g3.h5'
    }
  Nfiles = length(FileList);

  fprintf('Reshaping sst files:\n');

  for ifile = 1:Nfiles
    Fname = FileList{ifile};
 
    InFname = sprintf('%s/%s', SstPath, Fname);
    OutFname = sprintf('SST/%s', Fname);


    % Copy the input file datasets to the output, except for
    % SEATF (SST values) which needs to be reshaped before
    % writing out.
    fprintf(' Reading: %s\n', InFname);

    SEATF = h5read(InFname, '/SEATF');
    [ N1, N2 ] = size(SEATF);
    SEATF = reshape(SEATF, [ N2 N1 ]);

    DAY     = h5read(InFname, '/day'    );
    DX      = h5read(InFname, '/dx'     );
    DY      = h5read(InFname, '/dy'     );
    HOUR    = h5read(InFname, '/hour'   );
    MONTH   = h5read(InFname, '/month'  );
    NX      = h5read(InFname, '/nx'     );
    NY      = h5read(InFname, '/ny'     );
    POLELAT = h5read(InFname, '/polelat');
    POLELON = h5read(InFname, '/polelon');
    SW_LAT  = h5read(InFname, '/sw_lat' );
    SW_LON  = h5read(InFname, '/sw_lon' );
    YEAR    = h5read(InFname, '/year'   );
 

    fprintf(' Writing: %s\n', OutFname);
    if (exist(OutFname, 'file') == 2)
      delete(OutFname);
    end

    h5create(OutFname, '/SEATF',   size(SEATF)  );  
    h5create(OutFname, '/day',     size(DAY)    );
    h5create(OutFname, '/dx',      size(DX)     );
    h5create(OutFname, '/dy',      size(DY)     );
    h5create(OutFname, '/hour',    size(HOUR)   );
    h5create(OutFname, '/month',   size(MONTH)  );
    h5create(OutFname, '/nx',      size(NX)     );
    h5create(OutFname, '/ny',      size(NY)     );
    h5create(OutFname, '/polelat', size(POLELAT));
    h5create(OutFname, '/polelon', size(POLELON));
    h5create(OutFname, '/sw_lat',  size(SW_LAT) );
    h5create(OutFname, '/sw_lon',  size(SW_LON) );
    h5create(OutFname, '/year',    size(YEAR)   );

    h5write(OutFname, '/SEATF',   SEATF  );  
    h5write(OutFname, '/day',     DAY    );
    h5write(OutFname, '/dx',      DX     );
    h5write(OutFname, '/dy',      DY     );
    h5write(OutFname, '/hour',    HOUR   );
    h5write(OutFname, '/month',   MONTH  );
    h5write(OutFname, '/nx',      NX     );
    h5write(OutFname, '/ny',      NY     );
    h5write(OutFname, '/polelat', POLELAT);
    h5write(OutFname, '/polelon', POLELON);
    h5write(OutFname, '/sw_lat',  SW_LAT );
    h5write(OutFname, '/sw_lon',  SW_LON );
    h5write(OutFname, '/year',    YEAR   );
 
    fprintf('\n');
  end

end

