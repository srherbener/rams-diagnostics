function [ ] = TileFigs( FigList, TileDims, OutFile )
%
% TileFigs function that will take a list of MATLAB figure (.fig) files, generate
%          the images for these files and tile them according to given dimensions.
%
% Arguments:
%
%   FigList - cell array, list of the file names (with paths) of the
%             MATLAB figure files.
%
%   TileDims - vector, [ nrows, ncols ], specifying the tile organization
%              The figures given in FigFiles will fill out rows first, then
%              columns. Eg, if you give six files in FigFiles,
%              { 'f1' 'f2' 'f3' 'f4' 'f5' 'f6' } and spec TileDims as
%              [ 3 2 ], then the tiling of the images will go:
%
%                   f1 f2
%                   f3 f4
%                   f5 f6
%
%   OutFile - path to output file, the extension will determine what type of
%             image will be saved.
%

  % Use a built in format: JPEG, PNG or BMP for now
  Ofmt = 'NONE';
  Oext = 'NONE';
  if (~isempty(regexp(OutFile, '\.jpg')))
    Ofmt = '-djpeg';
    Oext = 'jpg';
  end
  if (~isempty(regexp(OutFile, '\.png')))
    Ofmt = '-dpng';
    Oext = 'png';
  end
  if (~isempty(regexp(OutFile, '\.bmp')))
    Ofmt = '-dbmp';
    Oext = 'bmp';
  end

  if (strcmp(Ofmt, 'NONE'))
    fprintf('TileArgs: ERROR: unrecognized file extension (image format) on OutFile: %s\n', OutFile);
    return;
  end

  % For raster images:
  %
  %   75 dpi -->  600 x  450 pixels
  %  150 dpi --> 1200 x  900 pixels
  %  300 dpi --> 2400 x 1800 pixels
  %  600 dpi --> 4800 x 3600 pixels
  %
  % Want the final pixel resolution to be 300dpi
  OutRes = 150;
  R = OutRes / 75;
  OutWidth = R * 600;
  OutHeight = R * 450;

  % Get a temporary file name that is unique so that multiple processes can
  % be running at the same time.
  TempFile = sprintf('%s.%s', tempname, Oext);
  
  % Read in images
  Nimg = length(FigList);
  
  Images = zeros(1,Nimg);
  for i = 1:Nimg
      Images(i) = openfig(FigList{i});
  end
  
  % Figure out how to size tiles.
  Dmax = max(TileDims);
  Ires =  OutRes / Dmax;
  Iwidth = OutWidth / Dmax;
  Iheight = OutHeight / Dmax;
  
  % Figure out the offsets for each tile within the output range
  k = 1;
  for j = TileDims(1):-1:1
      for i = 1:TileDims(2)
          Xtile(k) = (i-1) * Iwidth;
          Ytile(k) = (j-1) * Iheight;
          k = k + 1;
      end
  end
  
  % Do the tiling
  % The output needs to have the 
  Fig = figure('unit', 'pixel', 'position', [ 0 0 OutWidth OutHeight ], 'menubar', 'none');
  for i = 1:Nimg
      % Write out the image into a like format file as the output, then
      % read it in and place it in the output figure
      Ores = sprintf('-r%d', Ires);
      print(Images(i), Ofmt, Ores, TempFile);
      Img = imread(TempFile);
      AX = axes('parent', Fig, 'unit', 'pixel', 'position', [ Xtile(i) Ytile(i) Iwidth Iheight ]);
      imagesc(Img);
      axis(AX, 'off');
  end
  
  % Output the new image
  Ores = sprintf('-r%d', OutRes);
  print(Fig, Ofmt, Ores, OutFile);
  
  % Clean up
  for i = 1:Nimg
      close(Images(i));
  end
  close(Fig);

  delete(TempFile);
end

