function [ EofData ] = Xyzt2EofArray( Data4D, x1, x2, y1, y2, z1, z2, t1, t2)
%Azavg2EofArray Translate 4D data array to an array suitable for EOF
%analysis
%
%  This routine will take the 4D data, in the form (x,y,z,t), and
%  reformat it into an array that the EOF analysis can consume, which is in
%  the form (t, A) where A is a 1D representation of the 3D (x,y,z) field
%  created by the azavg routine. The format of the output array (EofData)
%  is:
%       rows --> time
%       columns --> 3D field
%           first Nx entries are  x = x1:x2, y = y1,   z = z1
%           second Nx entries are x = x1:x2, y = y1+1, z = z1
%           ...
%           Ny-th Nx entries are  x = x1:x2, y = y2, z = z1+1
%           next Nx entries are   x = x1:x2, y = y1, z = z1+1
%           ...


Nx = (x2-x1) + 1;
Ny = (y2-y1) + 1;
Nz = (z2-z1) + 1;
Nt = (t2-t1) + 1;

NA = Nx*Ny*Nz;
EofData = zeros(Nt, NA);

for it = t1:t2
    % draw out the (x,y,z) data at this time step into a linear (1D) array
    % only use the selected values for x,y,z however.
    iA = 1;
    A = zeros(1,NA);
    for iz = z1:z2
        for iy = y1:y2
            for ix = x1:x2
                A(iA) = Data4D(ix,iy,iz,it);
                iA = iA + 1;
            end
        end
    end
    EofData(it,:) = A;
end

end

