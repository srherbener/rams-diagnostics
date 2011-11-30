package DiagUtils;

use strict;
use File::Glob ':glob';
use File::Temp qw/ :mktemp  /;

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
  my %Plots;

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
      $Vars{$f[1]}{T_SUFFIX} = $f[3];
      }
    elsif ($f[0] eq "Azavg:")
      {
      $AzavgDiags{$f[1]}{DIM}    = $f[2];
      $AzavgDiags{$f[1]}{NBANDS} = $f[3];
      $AzavgDiags{$f[1]}{WMIN}   = $f[4];
      $AzavgDiags{$f[1]}{WMAX}   = $f[5];
      }
    elsif ($f[0] eq "Diag:")
      {
      $Diags{$f[1]}{$f[2]} = [ @f[3..$#f] ];
      }
    elsif ($f[0] eq "PlotExp:")
      {
      $Plots{EXP}      = $f[1];
      $Plots{STIME}    = $f[2];
      $Plots{TINC}     = $f[3];
      $Plots{PTSTART}  = $f[4];
      $Plots{PTEND}    = $f[5];
      }
    elsif ($f[0] eq "PlotTimes:")
      {
      $Plots{TIMES}  = [ @f[1..$#f] ];
      }
    elsif ($f[0] eq "PlotTypes:")
      {
      $Plots{TYPES}  = [ @f[1..$#f] ];
      }
    elsif ($f[0] eq "AzPlotVar:")
      {
      $Plots{"AZ"}{VARS}{$f[1]}{PTYPE} = $f[2];
      $Plots{"AZ"}{VARS}{$f[1]}{PYMIN} = $f[3];
      $Plots{"AZ"}{VARS}{$f[1]}{PYMAX} = $f[4];
      }
    elsif ($f[0] eq "HsPlotVar:")
      {
      # set a dummy value just to get the variable names into the keys
      $Plots{"HS"}{VARS}{$f[1]} = 1;
      }
    elsif ($f[0] eq "TsPlotVar:")
      {
      $Plots{"TS"}{VARS}{$f[1]}{TAVGLEN} = $f[2];
      $Plots{"TS"}{VARS}{$f[1]}{PYMIN}    = $f[3];
      $Plots{"TS"}{VARS}{$f[1]}{PYMAX}    = $f[4];
      }
    elsif ($f[0] eq "TsMplot:")
      {
      $Plots{"TSM"}{$f[1]}{$f[2]}{XLGD}  = $f[3];
      $Plots{"TSM"}{$f[1]}{$f[2]}{YLGD}  = $f[4];
      $Plots{"TSM"}{$f[1]}{$f[2]}{CASES} = [ @f[5..$#f] ];
      }
    }
  close(CONFIG);

  return(\%Cases, \%TimeDirs, \%Vars, \%AzavgDiags, \%Diags, \%Plots);
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
  my $FilePrefix;

  # append the new file name to the list
  if (($Var eq "cint_liq") || ($Var eq "cint_liq_sc") || ($Var eq "cint_liq_wr"))
    {
    $FilePrefix = "GRADS/" . $Tdir . "/" . $Var . "_" . $Case . "*";
    }
  else
    {
    $FilePrefix = $Case . "/GRADS/" . $Tdir . "/" . $Var . "-*";
    }

  ($InFile) = &FindGradsControlFile($FilePrefix);
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
  my ($FilePrefix) = @_;

  my $Gfile;
  my @f;

  # Assume that we will match (with bsd_glob) just the control file --> take
  # the first entry in the array that is returned from bsd_glob. 
  # Also assume that the GRADS control file we are looking for is one of
  # the "single variable" files where its name is of the form:
  #
  #    <variable>-<date_string>-<grid_number>.ctl
  $Gfile = $FilePrefix . ".ctl";
  @f = bsd_glob($Gfile);

  $Gfile = $f[0];

  return ($Gfile);
  }

#####################################################################################
# FixGradsControlFile()
#
# This routine will "fix" the given GRADS control file by:
#   1. replacing the path in the DSET command with just the file name so that
#      the control/data file pair will be portable
#   2. replacing the revu variable name with the given variable name
#
sub FixGradsControlFile
  {
  my ($GradsControlFile, $GradsDir, $Var) = @_; 

  my $BackupFile;
  my @SysArgs;

  my @f;
  my $NextOne;

  print "Fixing GRADS control file: $GradsControlFile\n";

  $BackupFile = $GradsControlFile . ".bak";
  @SysArgs = ("mv", $GradsControlFile, $BackupFile);
  system(@SysArgs);

  open(BACKUP, $BackupFile) or die "Cannot open $BackupFile for reading: $!";
  open(G_CTRL, ">$GradsControlFile") or die "Cannot open $GradsControlFile for writing: $!";

  # assume that there is only one variable, so look for the vars keyword and make the
  # change on the next line (following the "vars" line)
  $NextOne = 0;
  while (<BACKUP>)
    {   
    if (/^[Dd][Ss][Ee][Tt]/)
      {
      s/$GradsDir\///;
      }
    elsif (/^[Vv][Aa][Rr][Ss]/)
      {
      $NextOne = 1;
      }
    elsif ($NextOne == 1)
      {
      @f = split(' ');
      s/\b$f[0]\b/$Var/;
      $NextOne = 0;
      }

    print G_CTRL;
    }   
  close(BACKUP);
  close(G_CTRL);

  return;
  }

#####################################################################################
# PlotGradsVslice()
#
# This routine will run grads and create a plot, vertical slice through the data,
# using the specs given in the arguments.
#

sub PlotGradsVslice
  {
  my ($Gvar, $Gfile, $TimeStep, $Z1, $Z2, $Ccols, $Clevs, $Ptitle, $Xtitle, $Ytitle, $Pfile) = @_;

  my $GscriptName;
  my $GscriptFh;

  my $Var15; # GRADS truncates the variable names to 15 characters

  $Var15 = substr($Gvar, 0, 15); # copy the first 15 characters

  # mkstemp opens the file
  ($GscriptFh, $GscriptName) = mkstemp( "/tmp/gcmds.XXXXXXXX" );
  print $GscriptFh "export GASCRP=\"\$HOME/grads\"\n";
  print $GscriptFh "grads -l -b <<EOF\n";
  print $GscriptFh "reinit\n";
  print $GscriptFh "clear\n";
  print $GscriptFh "open $Gfile\n";
  print $GscriptFh "set lev $Z1 $Z2\n";
  print $GscriptFh "set t $TimeStep\n";
  print $GscriptFh "set gxout shaded\n";
  print $GscriptFh "set clevs $Clevs\n";
  print $GscriptFh "set ccols $Ccols\n";
  print $GscriptFh "set xlab %.1f\n";
  print $GscriptFh "set grads off\n";
  print $GscriptFh "set parea 1.5 10.5 2 8\n";
  print $GscriptFh "d $Var15\n";
  print $GscriptFh "draw title $Ptitle\n";
  print $GscriptFh "draw xlab $Xtitle\n";
  print $GscriptFh "draw ylab $Ytitle\n";
  print $GscriptFh "cbarn 1.0 0\n";
  print $GscriptFh "printim $Pfile white\n";
  print $GscriptFh "EOF\n";
  close($GscriptFh);

  system ("chmod", "+x", $GscriptName);
  system ("bash" , "-c", $GscriptName);
 
  unlink $GscriptName or warn "PlotGradsVslice: could not unlink $GscriptName: $!";

  return;
  }

#####################################################################################
# PlotGradsHslice()
#
# This routine will run grads and create a plot, horizontal slice through the data,
# using the specs given in the arguments.
#

sub PlotGradsHslice
  {
  my ($Gvar, $Gfile, $TimeStep, $Ccols, $Clevs, $Ptitle, $Xtitle, $Ytitle, $Pfile) = @_;

  my $GscriptName;
  my $GscriptFh;

  my $Var15; # GRADS truncates the variable names to 15 characters

  $Var15 = substr($Gvar, 0, 15); # copy the first 15 characters

  # mkstemp opens the file
  ($GscriptFh, $GscriptName) = mkstemp( "/tmp/gcmds.XXXXXXXX" );
  print $GscriptFh "export GASCRP=\"\$HOME/grads\"\n";
  print $GscriptFh "grads -l -b <<EOF\n";
  print $GscriptFh "reinit\n";
  print $GscriptFh "clear\n";
  print $GscriptFh "open $Gfile\n";
  print $GscriptFh "set t $TimeStep\n";
  print $GscriptFh "set gxout shaded\n";
  print $GscriptFh "set clevs $Clevs\n";
  print $GscriptFh "set ccols $Ccols\n";
  print $GscriptFh "set grads off\n";
  print $GscriptFh "d $Var15\n";
  print $GscriptFh "draw title $Ptitle\n";
  print $GscriptFh "draw xlab $Xtitle\n";
  print $GscriptFh "draw ylab $Ytitle\n";
  print $GscriptFh "cbarn 1.0 1\n";
  print $GscriptFh "printim $Pfile white\n";
  print $GscriptFh "EOF\n";
  close($GscriptFh);

  system ("chmod", "+x", $GscriptName);
  system ("bash" , "-c", $GscriptName);
 
  unlink $GscriptName or warn "PlotGradsHslice: could not unlink $GscriptName: $!";

  return;
  }

#####################################################################################
# PlotGradsTseries()
#
# This routine will run grads and create a plot, horizontal slice through the data,
# using the specs given in the arguments.
#

sub PlotGradsTseries
  {
  my ($Gvar, $Gfiles, $CaseList, $TstepStart, $TstepEnd, $TavgLen, $Ymin, $Ymax,
      $Ptitle, $Xtitle, $Ytitle, $Xlgd, $Ylgd, $Pfile) = @_;

  my $GscriptName;
  my $GscriptFh;
  my $T1;
  my $T2;
  my $Var15; # GRADS truncates the variable names to 15 characters
  my $i;
  my $ip1;
  my $CaseString;

  $T1 = $TstepStart + $TavgLen;
  $T2 = $TstepEnd - $TavgLen;

  $Var15 = substr($Gvar, 0, 15); # copy the first 15 characters
  $CaseString = '"' . $$CaseList[0] . '"';

  # mkstemp opens the file
  ($GscriptFh, $GscriptName) = mkstemp( "/tmp/gcmds.XXXXXXXX" );
  print $GscriptFh "export GASCRP=\"\$HOME/grads\"\n";
  print $GscriptFh "grads -l -b <<EOF\n";
  print $GscriptFh "reinit\n";
  print $GscriptFh "clear\n";
  print $GscriptFh "open $$Gfiles[0]\n";
  print $GscriptFh "set t $T1 $T2\n";
  print $GscriptFh "set grads off\n";
  print $GscriptFh "set vrange $Ymin $Ymax\n";
  print $GscriptFh "d tloop(ave($Var15,t-$TavgLen,t+$TavgLen))\n";
  print $GscriptFh "draw title $Ptitle\n";
  print $GscriptFh "draw xlab $Xtitle\n";
  print $GscriptFh "draw ylab $Ytitle\n";
  foreach $i (1..$#$Gfiles)
    {
    $ip1 = $i + 1;
    print $GscriptFh "open $$Gfiles[$i]\n";
    print $GscriptFh "set dfile $ip1\n";
    print $GscriptFh "d tloop(ave($Var15,t-$TavgLen,t+$TavgLen))\n";
    $CaseString = $CaseString . ' "' . $$CaseList[$i] . '"';
    }
  if ($#$Gfiles > 0)
    {
    print $GscriptFh "cbarline -x $Xlgd -y $Ylgd -t $CaseString\n";
    }
  print $GscriptFh "printim $Pfile white\n";
  print $GscriptFh "EOF\n";
  close($GscriptFh);

  system ("chmod", "+x", $GscriptName);
  system ("bash" , "-c", $GscriptName);
 
  unlink $GscriptName or warn "PlotGradsTslice: could not unlink $GscriptName: $!";

  return;
  }

#####################################################################################
# PlotGradsHovmol()
#
# This routine will run grads and create a plot, vertical slice through the data,
# using the specs given in the arguments.
#

sub PlotGradsHovmol
  {
  my ($Gvar, $Gfile, $TstepStart, $TstepEnd, $Ccols, $Clevs,
      $Ptitle, $Xtitle, $Ytitle, $Pfile) = @_;

  my $GscriptName;
  my $GscriptFh;

  my $Var15; # GRADS truncates the variable names to 15 characters

  $Var15 = substr($Gvar, 0, 15); # copy the first 15 characters

  # mkstemp opens the file
  ($GscriptFh, $GscriptName) = mkstemp( "/tmp/gcmds.XXXXXXXX" );
  print $GscriptFh "export GASCRP=\"\$HOME/grads\"\n";
  print $GscriptFh "grads -l -b <<EOF\n";
  print $GscriptFh "reinit\n";
  print $GscriptFh "clear\n";
  print $GscriptFh "open $Gfile\n";
  print $GscriptFh "set t $TstepStart $TstepEnd\n";
  print $GscriptFh "set gxout shaded\n";
  print $GscriptFh "set clevs $Clevs\n";
  print $GscriptFh "set ccols $Ccols\n";
  print $GscriptFh "set xlab %.1f\n";
  print $GscriptFh "set grads off\n";
  print $GscriptFh "set parea 1.5 10.5 2 8\n";
  print $GscriptFh "d $Var15\n";
  print $GscriptFh "draw title $Ptitle\n";
  print $GscriptFh "draw xlab $Xtitle\n";
  print $GscriptFh "draw ylab $Ytitle\n";
  print $GscriptFh "cbarn 1.0 0\n";
  print $GscriptFh "printim $Pfile white\n";
  print $GscriptFh "EOF\n";
  close($GscriptFh);

  system ("chmod", "+x", $GscriptName);
  system ("bash" , "-c", $GscriptName);
 
  unlink $GscriptName or warn "PlotGradsHovmol: could not unlink $GscriptName: $!";

  return;
  }

#####################################################################################
# GetGradsStats()
#
# This routine will run grads and load stats for a given variable in a given
# GRADS file into a hash that gets returned to the caller.
#

sub GetGradsStats
  {
  my ($Gfile, $Gvar) = @_;

  my %Gstats;

  my $GscriptName;
  my $GscriptFh;
  my $Cmd;
  my @f;
  my $Time;
  my $Zlev;
  my $Min;
  my $Max;
  my $Mean;
  my $StdDev;

  my $Var15; # GRADS truncates the variable names to 15 characters

  $Var15 = substr($Gvar, 0, 15); # copy the first 15 characters
  undef(%Gstats);

  # mkstemp opens the file
  ($GscriptFh, $GscriptName) = mkstemp( "/tmp/gcmds.XXXXXXXX" );
  print $GscriptFh "export GASCRP=\"\$HOME/grads\"\n";
  print $GscriptFh "grads -l -b <<EOF\n";
  print $GscriptFh "reinit\n";
  print $GscriptFh "clear\n";
  print $GscriptFh "dumpStats $Gfile $Var15\n";
  print $GscriptFh "EOF\n";
  close($GscriptFh);

  system ("chmod", "+x", $GscriptName);

  # run dumpStats (GRADS script) and collect output in %Gstats
  # look for output like this:
  #   Stats:
  #     t: 1
  #     z: 1
  #     Min: 0
  #     Max: 1.41709
  #     Mean: 0.161222
  #     StdDev: 0.154379

  $Cmd = "bash -c $GscriptName";
  open (STATS, "$Cmd |") or die "GetGradsStats: cannot run command: $Cmd: $!";
  while (<STATS>)
    {
    @f = split(' ');
    if ($f[0] eq "t:")
      {
      $Time = $f[1];
      }
    elsif ($f[0] eq "z:")
      {
      $Zlev = $f[1];
      }
    elsif ($f[0] eq "Min:")
      {
      $Min = $f[1];
      }
    elsif ($f[0] eq "Max:")
      {
      $Max = $f[1];
      }
    elsif ($f[0] eq "Mean:")
      {
      $Mean = $f[1];
      }
    elsif ($f[0] eq "StdDev:")
      {
      $StdDev = $f[1];
      $Gstats{$Time}{$Zlev}{MIN} = $Min;
      $Gstats{$Time}{$Zlev}{MAX} = $Max;
      $Gstats{$Time}{$Zlev}{MEAN} = $Mean;
      $Gstats{$Time}{$Zlev}{STD_DEV} = $StdDev;
      }
    }
  close(STATS);
 
  unlink $GscriptName or warn "GetGradsStats: could not unlink $GscriptName: $!";

  return (\%Gstats);
  }

1;
