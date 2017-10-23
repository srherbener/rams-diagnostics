; IDL procedure to generate plots
;
; This routine depend on compiling the read_grads and plot_utils routines before running.
;

; put this in so we can index arrays beyond 32,768 entries
Compile_Opt DEFINT32

;**********************************************************
; read_revu_hdf5_var
;
; Args
;  1. hdf5 file
;  2. grid number
;  3. variable name
;  4, output variable data
;

pro read_revu_hdf5_var, file, grid, vname, out_var

  dset_path = string(format='("/g",i0,"/",a)', grid, vname)

  ; open and read in the variable data
  fileid = H5F_OPEN(file)
  dsetid = H5D_OPEN(fileid,dset_path)
  out_var = H5D_READ(dsetid)

  ; out_var has data from hdf5 file, but dimensions are reversed
  out_var = transpose(out_var)

  H5D_CLOSE, dsetid
  H5F_CLOSE, fileid

end
