package DiagUtils;

use strict;
use File::Glob ':glob';

###############################################################################
# ReadConfigFile()
#
# This routine will read in info from the diagnostics configuration file.
# Hashes will be loaded up with the information and passed back to the caller.
#

sub ReadConfigFile
  {
  my ($ConfigFile) = @_;

  my %Cases;
  my %TimeDirs;
  my %Vars;
  my %AzavgDiags;
  my %Diags;

  my @f;

  undef(%Cases);
  undef(%TimeDirs);
  undef(%Vars);
  undef(%AzavgDiags);

  open(CONFIG, "$ConfigFile") or die "Cannot open '$ConfigFile' for reading: $!";
  while (<CONFIG>)
    {
    @f = split(' ');
  
    if ($f[0] eq "Case:")
      {
      $Cases{$f[1]} = 1;
      }
    elsif ($f[0] eq "TimeDir:")
      {
      $TimeDirs{$f[1]} = $f[2];
      }
    elsif ($f[0] eq "Var:")
      {
      $Vars{$f[1]}{REVU_VAR} = $f[2];
      $Vars{$f[1]}{GRADS_VAR} = $f[3];
      $Vars{$f[1]}{T_SUFFIX} = $f[4];
      }
    elsif ($f[0] eq "Azavg:")
      {
      $AzavgDiags{$f[1]}{NBANDS} = $f[2];
      $AzavgDiags{$f[1]}{WMIN} = $f[3];
      $AzavgDiags{$f[1]}{WMAX} = $f[4];
      }
    elsif ($f[0] eq "Diag:")
      {
      $Diags{$f[1]}{$f[2]} = [ @f[3..$#f] ];
      }
    }
  close(CONFIG);

  return(\%Cases, \%TimeDirs, \%Vars, \%AzavgDiags, \%Diags);
  }

#######################################################################
# AppendGradsVarFile()
#
# This routine will take the case, time dir, and diagnostic name
# and create the file lists (input and output) needed for azavg.
#
sub AppendGradsVarFile
  {
  my ($InFiles, $Case, $Tdir, $Var) = @_;

  my $InFileList;

  my $InFile;
  my @f;

  # append the new file name to the list
  ($InFile) = &FindGradsControlFile($Case, $Tdir, $Var);
  $InFileList = $InFiles;
  if ($InFileList eq "")
    {
    $InFileList = $InFile;
    }
  else
    {
    $InFileList = $InFileList . ":" . $InFile;
    }
 
  return $InFileList;
  }

###############################################################################
# FindGradsControlFile()
#
# This routine will take in a case, time dir and var name and return the
# path to the corresponding GRADS control file
#

sub FindGradsControlFile
  {
  my ($Case, $Tdir, $Var) = @_;

  my $Gfile;
  my @f;

  # Assume that we will match (with bsd_glob) just the control file --> take
  # the first entry in the array that is returned from bsd_glob. 
  # Also assume that the GRADS control file we are looking for is one of
  # the "single variable" files where its name is of the form:
  #
  #    <variable>-<date_string>-<grid_number>.ctl
  $Gfile = $Case . "/GRADS/" . $Tdir . "/" . $Var . "-*.ctl";
  @f = bsd_glob($Gfile);

  $Gfile = $f[0];

  return ($Gfile);
  }

1;
