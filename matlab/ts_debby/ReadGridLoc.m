function [ Gloc ] = ReadGridLoc(Fname, Glabel)

  InFid = fopen(Fname);

  % Read entire file into cell array
  % Each line of the file is: Label lat1 lat2 lon1 lon2
  % Format for textscan reflects that for one line
  Gdata = textscan(InFid, '%s %f %f %f %f');

  % Each entry in Gdata is a column from the file
  Labels = Gdata{1};
  Locs = [ Gdata{2} Gdata{3} Gdata{4} Gdata{5} ];

  Gloc = [ 0 0 0 0 ];
  for i = 1:length(Labels)
    if (strcmp(Labels{i}, Glabel))
      Gloc = squeeze(Locs(i,:));
    end
  end

  fclose(InFid);
end
