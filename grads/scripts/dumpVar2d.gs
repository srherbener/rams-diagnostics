************************************************************************
* dump out ascii for a 2D variable (X-Z vertical slices)
*   each line is a time step
*   the entries on a line contain the variable value at each x,y location
*   the first entry is the lower left corner, the last entry the upper right
*
* varNum arg:
*   1 --> var file
*   2 --> var name
*   3 --> number of time points
*   4 --> number of x points
*   5 --> number of y points

function main(args)

  Gfile  = subwrd(args, 1)
  Gvar   = subwrd(args, 2)
  Tstart = subwrd(args, 3)
  Tend   = subwrd(args, 4)
  Xstart = subwrd(args, 5)
  Xend   = subwrd(args, 6)
  Zstart = subwrd(args, 7)
  Zend   = subwrd(args, 8)
  Ofile  = subwrd(args, 9)

  'reinit'
  'clear'
  'open 'Gfile

  it = Tstart
  while (it <= Tend)
    'set t 'it
    OutLine = ''
    iz = Zstart
    while (iz <= Zend)
      'set z 'iz
      ix = Xstart
      while (ix <= Xend)
        'set x 'ix

        'd 'Gvar
        Value = subwrd(result,4)
        if (Value = -9.99e+08)
          Value = 'NaN'
        endif
        OutLine = OutLine' 'Value

        ix = ix + 1
      endwhile
      iz = iz + 1
    endwhile
    RetCode = write(Ofile, OutLine)
    it = it + 1
  endwhile

  RetCode = close(Ofile)

return
