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

  my %Config;

  my @f;

  undef(%Config);

  open(CONFIG, "$ConfigFile") or die "Cannot open '$ConfigFile' for reading: $!";
  while (<CONFIG>)
    {
    @f = split(' ');
  
    if ($f[0] eq "Case:")
      {
      $Config{CASES}{$f[1]} = 1;
      }
    elsif ($f[0] eq "RamsDir:")
      {
      $Config{RAMS_DIR} = $f[1];
      }
    elsif ($f[0] eq "RevuDir:")
      {
      $Config{REVU_DIR} = $f[1];
      }
    elsif ($f[0] eq "DiagDir:")
      {
      $Config{DIAG_DIR} = $f[1];
      }
    elsif ($f[0] eq "AzavgDir:")
      {
      $Config{AZAVG_DIR} = $f[1];
      }
    elsif ($f[0] eq "TsavgDir:")
      {
      $Config{TSAVG_DIR} = $f[1];
      }
    elsif ($f[0] eq "FilterDir:")
      {
      $Config{FILTER_DIR} = $f[1];
      }
    elsif ($f[0] eq "PlotDir:")
      {
      $Config{PLOT_DIR} = $f[1];
      }
    elsif ($f[0] eq "EofDir:")
      {
      $Config{EOF_DIR} = $f[1];
      }
    elsif ($f[0] eq "TimeDir:")
      {
      $Config{TIM_DIRS}{$f[1]} = $f[2];
      }
    elsif ($f[0] eq "Var:")
      {
      $Config{VARS}{$f[1]}{REVU_VAR} = $f[2];
      $Config{VARS}{$f[1]}{T_SUFFIX} = $f[3];
      }
    elsif ($f[0] eq "Azavg:")
      {
      $Config{AZAVG_DIAGS}{$f[1]}{VAR}       = $f[2];
      $Config{AZAVG_DIAGS}{$f[1]}{DIM}       = $f[3];
      $Config{AZAVG_DIAGS}{$f[1]}{NBANDS}    = $f[4];
      $Config{AZAVG_DIAGS}{$f[1]}{FILTER}    = $f[5];
      }
    elsif ($f[0] eq "Tsavg:")
      {
      $Config{TSAVG_DIAGS}{$f[1]}{DSPEC} = $f[2];
      $Config{TSAVG_DIAGS}{$f[1]}{FILTER} = $f[3];
      }
    elsif ($f[0] eq "Op:")
      {
      $Config{OPS}{$f[1]}{VAR}   = $f[2];
      $Config{OPS}{$f[1]}{CASE1} = $f[3];
      $Config{OPS}{$f[1]}{CASE2} = $f[4];
      $Config{OPS}{$f[1]}{OP}    = $f[5];
      }
    elsif ($f[0] eq "Filter:")
      {
      $Config{FILTERS}{$f[1]}{MODEL} = $f[2];
      $Config{FILTERS}{$f[1]}{SPECS} = [ @f[3..$#f] ];
      }
    elsif ($f[0] eq "Moment:")
      {
      $Config{MOMENTS}{$f[1]}{FILTER} = $f[2];
      $Config{MOMENTS}{$f[1]}{SPECS}  = [ @f[3..$#f] ];
      }
    elsif ($f[0] eq "Diag:")
      {
      $Config{DIAGS}{$f[1]}{$f[2]} = [ @f[3..$#f] ];
      }
    }
  close(CONFIG);

  return(\%Config);
  }

1;
