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

  % Use a built in format: PNG, JPEG
  Ofmt = '-djpeg';
  Oext = 'jpg';
  
  % Want the final pixel size to be 1200 x 900 (150dpi)
  OutWidth = 1200;
  OutHeight = 900;
  OutRes = 150;

  % Get a temporary file name that is unique so that multiple processes can
  % be running at the same time.
  TempFile = sprintf('%s.%s', tempname, Oext);
  
  % Read in images
  Nimg = length(FigList);
  
  Images = zeros(1,Nimg);
  for i = 1:Nimg
      Images(i) = openfig(FigList{i});
  end
  
  % Figure out how to size tiles. Really only supporting dimensions sizes
  % of 1, 2 or 3 for now.
  %
  % For JPEG and PNG:
  %   75 dpi --> 600 x 451 pixels
  %  150 dpi --> 1200 x 900 pixels
  %
  % Want the resulting format to end up at 1200 x 900 (150dpi) so 75dpi
  % is good for 2x2, 50dpi for 3x3, etc.
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

