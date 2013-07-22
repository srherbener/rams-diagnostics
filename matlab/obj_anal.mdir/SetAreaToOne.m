function [ NewY ] = SetAreaToOne( X, Y )
%SetAreaToOne Set the aree under the curve to one.
%   This function takes the x and y values and adjusts the y values so that
%   the aree under the curve becomes equal to one.

LenX = length(X);
LenY = length(Y);
if (LenX == LenY)
    % lengths must be equal
    TotalArea = 0;
    NewY = zeros(1,LenX);
    for i = 2:LenX;
        TotalArea = TotalArea + (Y(i-1) * (X(i)-X(i-1)));
    end
    NewY = Y / TotalArea;
else
    fprintf('ERROR: SetAreaToOne: length of X and Y vectors must be equal');
    NewY = nan(1,LenY);
end

end

