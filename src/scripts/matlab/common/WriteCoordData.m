function [ ] = WriteCoordData(Coords, OutFile)
% WriteCooordData function to write out coordinate values into a text file

  file_id = fopen(OutFile, 'w');
  for i = 1:length(Coords)
    fprintf(file_id, ' %.3f', Coords(i));
  end 
  fprintf(file_id, '\n');
  fclose(file_id);
end

