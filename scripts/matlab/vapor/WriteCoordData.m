function WriteCoordData(data, OutFile)
% WriteCooordData function to write out coordinate values into a text file

    fID = fopen(OutFile,'w');

    for n=1:length(data)
        fprintf( fID, ' %.3f', data(n) );
    end
    fprintf( fID, '\n' );

    fclose( fID );

end
