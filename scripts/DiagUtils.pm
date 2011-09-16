package DiagUtils;

use strict;

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

1;
