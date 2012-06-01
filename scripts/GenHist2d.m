function [ Hists ] = GenHist2d( Var, Bins, Hmin, Hmax, Scale )
%GenHist2d create histograms of 2D field
%   This routine will create histograms of the series of 2D fields given
%   by Var. It is assumed that Var has 3 dimensions organized as (x,y,t).
%   The x,y portions are the 2D fields, and t is the "series" dimension.
%   One histogram is generated per series and passed back thorugh Hists
%   which is a 2D array organized as (hist, t), ie. each column is a
%   histogram.
%
%   The argument Bins defines the bins for the histogram in that
%   each value in Bins is an edge for that bin. Edge means that if the
%   value you are testing, x, satisfies Bins(i) <= x < Bins(i+1), then the
%   count is incremented for element i in the output histogram. The
%   histogram will have the same length as the vector Bins since the last
%   bin holds a count of the x's that equal Bins(n).
%
%   The arguments Hmin, Hmax define a range for selecting items out of Var
%   where only the elements in Var that fall inside the range Hmin, Hmax
%   will be counted for the output histogram.
%
%   The argument Scale defines how to scale the resulting histograms.
%      'FA' --> fractional area

[ Nx, Ny, Nt ] = size(Var);
Nbins = length(Bins);
Npts = Nx * Ny;

Hists = zeros(Nbins, Nt); % preallocate to help speed up things

for i = 1:Nt
    % Pull out the 2D field
    Field = squeeze(Var(:,:,i));
    
    % Apply the selection
    Hdata = Field(Field >= Hmin & Field <= Hmax);
    
    % Generate the histogram
    H = histc(Hdata, Bins);
    
    % Attach to the output array
    Hists(:,i) = H;
end

if (Scale == 'FA')
  Hists = Hists / Npts;
end

end

