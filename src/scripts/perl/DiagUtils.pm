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
      $Config{AZAVG_DIAGS}{$f[1]}{V_SPECS}     = $f[2];
      $Config{AZAVG_DIAGS}{$f[1]}{V_NAME}      = $f[3];
      $Config{AZAVG_DIAGS}{$f[1]}{NUM_RBANDS}  = $f[4];
      $Config{AZAVG_DIAGS}{$f[1]}{MAX_RADIUS}  = $f[5];
      $Config{AZAVG_DIAGS}{$f[1]}{F_PREFIX}    = $f[6];
      $Config{AZAVG_DIAGS}{$f[1]}{SC_PREFIX}   = $f[7];
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
      $Config{FILTERS}{$f[1]}{SPECS} = [ @f[2..$#f] ];
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

###############################################################################
# ReadRevuConfigFile()
#
# This routine will read in info from the REVU configuration file.
# Hashes will be loaded up with the information and passed back to the caller.
#

sub ReadRevuConfigFile
  {
  my ($ConfigFile) = @_;

  my %Config;

  my @f;

  undef(%Config);

  $Config{CFORMAT} = "new";

  open(CONFIG, "$ConfigFile") or die "Cannot open '$ConfigFile' for reading: $!";
  while (<CONFIG>)
    {
    @f = split(' ');
  
    if ($f[0] eq "AnPref:")
      {
      $Config{ANPREF} = $f[1];
      }
    elsif ($f[0] eq "RevPref:")
      {
      $Config{REVPREF} = $f[1];
      }
    elsif ($f[0] eq "AnaType:")
      {
      $Config{ANATYPE} = $f[1];
      }
    elsif ($f[0] eq "Igrid:")
      {
      $Config{IGRID} = $f[1];
      }
    elsif ($f[0] eq "Iztran:")
      {
      $Config{IZTRAN} = $f[1];
      }
    elsif ($f[0] eq "Case:")
      {
      $Config{CASES}{$f[1]} = 1;
      }
    elsif ($f[0] eq "Cformat:")
      {
      $Config{CFORMAT} = $f[1];
      }
    elsif ($f[0] eq "Var:")
      {
      $Config{VARS}{$f[1]}{REVU_VAR} = $f[2];
      $Config{VARS}{$f[1]}{XVAR}     = $f[3];
      $Config{VARS}{$f[1]}{YVAR}     = $f[4];
      $Config{VARS}{$f[1]}{ZVAR}     = $f[5];
      $Config{VARS}{$f[1]}{TVAR}     = $f[6];
      }
    }
  close(CONFIG);

  return(\%Config);
  }


###############################################################################
# ReadDiagConfigFile()
#
# This routine will read in info from the Diagnostics configuration file, and
# organize the information for the caller.
#

sub ReadDiagConfigFile
  {
  my ($ConfigFile) = @_;

  my %Config;

  my $Cstmts;
  my $i;
  my $j;
  my @f;

  # reformat the input file into a list of statements
  ($Cstmts) = &PrepConfigFile($ConfigFile);

  foreach $i (0 .. $#$Cstmts)
    {
    @f = @{ $$Cstmts[$i] };

    if ($f[0] eq "CaseSet:")
      {
      $Config{CASE_LIST} = [ @f[ 1 .. $#f ] ];
      }
    elsif ($f[0] eq "Azavg:")
      {
      $Config{AZAVG}{$f[1]}{NUM_RBANDS}  = $f[2];
      $Config{AZAVG}{$f[1]}{MAX_RADIUS}  = $f[3];
      $Config{AZAVG}{$f[1]}{DO_HIST}     = $f[4];
      $Config{AZAVG}{$f[1]}{SUB_SMOTION} = $f[5];
      $Config{AZAVG}{$f[1]}{IN_TYPE}     = $f[6];
      $Config{AZAVG}{$f[1]}{FILE_LIST}   = [ @f[7..$#f] ];
      }
    elsif ($f[0] eq "Tsavg:")
      {
      $Config{TSAVG}{$f[1]}{AVG_FUNC}    = $f[2];
      $Config{TSAVG}{$f[1]}{SUB_SMOTION} = $f[3];
      $Config{TSAVG}{$f[1]}{FILE_LIST}   = [ @f[4..$#f] ];
      }
    elsif ($f[0] eq "Filter:")
      {
      $Config{FILTER}{$f[1]}{COMB_SENSE}    = $f[2];
      $Config{FILTER}{$f[1]}{SPEC_LIST}   = [ @f[3..$#f] ];
      }
    elsif ($f[0] eq "Xsection:")
      {
      $Config{XSECTION}{$f[1]}{IN_TYPE}     = $f[2];
      $Config{XSECTION}{$f[1]}{LINE_SPEC}   = $f[3];
      $Config{XSECTION}{$f[1]}{FILE_LIST}   = [ @f[4..$#f] ];
      }
    elsif ($f[0] eq "Advect:")
      {
      $Config{ADVECT}{$f[1]}{ZMIN}        = $f[2];
      $Config{ADVECT}{$f[1]}{ZMAX}        = $f[3];
      $Config{ADVECT}{$f[1]}{FILE_LIST}   = [ @f[4..$#f] ];
      }
    elsif ($f[0] eq "VintTerms:")
      {
      $Config{VINT_TERMS}{$f[1]}{FILE_LIST}   = [ @f[2..$#f] ];
      }
    elsif ($f[0] eq "VtVr:")
      {
      $Config{VTVR}{$f[1]}{SUB_SMOTION} = $f[2];
      $Config{VTVR}{$f[1]}{REV_CONV}    = $f[3];
      $Config{VTVR}{$f[1]}{FILE_LIST}   = [ @f[4..$#f] ];
      }
    elsif ($f[0] eq "SubVort:")
      {
      $Config{SVORT}{$f[1]}{FILE_LIST}   = [ @f[2..$#f] ];
      }
  }

  return(\%Config);
  }

###############################################################################
# PrepConfigFile()
#
# This routine will read in a configuration file and reformat it into a list
# of statements.
#

sub PrepConfigFile
  {
  my ($ConfigFile) = @_;

  my @Cstmts;

  my @f;
  my $istmnt;
  my $itoken;
  my $i;


  $istmnt = 0;
  $itoken = 0;
  open(CONFIG, "$ConfigFile") or die "Cannot open $ConfigFile for reading: $!";
  while(<CONFIG>)
    {
    # strip off comments
    s/#.*//;

    # split into tokens and load these into the current statement
    # if we encounter the token 'End' then move to the next statement
    @f = split(' ');
    foreach $i (0..$#f)
      {
      if ($f[$i] eq "End")
        {
        $istmnt++;
        $itoken = 0;
        }
      else
        {
        $Cstmts[$istmnt][$itoken] = $f[$i];
        $itoken++;
        }
      }
    }
  close(CONFIG); 

  return (\@Cstmts);
  }

#############################################################################
# SubmitParallelJobs()
#
# This routine will take a list of commands (the jobs to be run), divide
# then up evenly among the given number of processes, and submit the jobs.
#

sub SubmitParallelJobs
  {
  my ($JobList, $Nprocs, $Dname) = @_;

  my $icmd;
  my @SysArgs;
  my $ilist;
  my @SubLists;

  my $ilog;
  my $LogFile;
  my $Pid;

  # Split up list into Nprocs sub-lists making the
  # sizes of these lists equal (some will have one
  # more item than others).
  foreach $icmd (0 .. $#{ $JobList })
    {
    $ilist = $icmd % $Nprocs;
    push @{ $SubLists[$ilist]}, [ @{ $$JobList[$icmd] } ];
    }

  print "Submitting diagnostic jobs on $Nprocs processes for $Dname\n";
  print "\n";

  # Run each sublist on a separate process.
  $ilog = 1;
  foreach $ilist (0 .. $#SubLists)
    {
    # Create numbered log files so each child process gets a unique
    # one. However, always number from 1 to num_cases so that a 
    # massive number of files doesn't collect after running this
    # command a bunch of times. The user will have to understand that
    # the log files get overwritten each time you run this command.

    $LogFile = sprintf("%s%d.log", $Dname, $ilog);

    $Pid = fork;

    # Child will run through the list of diagnostics and
    # quit. The parent will continue in the Case loop
    # (Pid == 0 for the parent).

    if ($Pid != 0)
      {
      # Child
  
      # Redirect STDOUT and STDERR into the log file
      open STDOUT, ">",  $LogFile;
      open STDERR, ">>", $LogFile;

      print "Running diagnostics on child process: $Pid\n";
      print "\n";

      foreach $icmd (0 .. $#{ $SubLists[$ilist] })
        {
        @SysArgs = @{ $SubLists[$ilist][$icmd] };
        print "Running: ", join(" ", @SysArgs), "\n";

        system(@SysArgs);
  
        print "\n";
        }

      exit 0;
      }

    # If made it to here, in the parent process
    $ilog++;
    }

  return;
  }

1;
