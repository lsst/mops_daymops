# jmyers i've pretty much forgotten how to write/run perl so i'm just
# copy-pasting from old PS mops code...

use Carp;
use POSIX;

our $MPC_FORMAT_SPEC = "%5.5s%-7.7s%1.1s%1.1s%1.1s%-17.17s%-12.12s%-12.12s         %-5.2f%1.1s      %-3.3s";
$StrZero = 0;


sub miti_read {
    # Read in MITI file. pair file.  Return list of ARRAYREFs.  If filename is "-", use
    # STDIN.
    my $filename = shift;
    my $line;
    my @items;

    open my $fh, $filename or die "can't open MITI file $filename";
    while (defined($line = <$fh>)) {
        if ($line !~ /^#/) {    # ignore comments
            push @items, [split /\s+/, $line];
        }
    }
    close $fh;
    return @items;
}
*miti_slurp = \&miti_read;





sub miti_parse {
    # Parse a string into MITI fields.  Return a hash containing the keys:
    #  ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME
    my $line = shift;
    my %result;
    my @items = split /\s+/, $line;     # separate tokens
    #print "line = $line, items = @items \n";
    @result{qw(ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME)} = @items;  # bulk assign
    return %result;
}




sub deg2str ($$$;$) {

  my($deg, $mode, $sig, $strsep) = @_;
  return turn2str($deg/(360), $mode, $sig, $strsep);
}




# Return the nearest integer (ie round)
sub nint ($) {
  my ($x) = @_;
  ($x<0.0) ? return(ceil($x-0.5)) : return(floor($x+0.5))
}





sub turn2str ($$$;$) {
  my($turn, $mode, $sig, $strsep) = @_;
  $mode = uc $mode;
  if (($mode ne 'H') && ($mode ne 'D')) {
    carp 'turn2str: $mode must equal \'H\' or \'D\'';
    return undef;
  }
  $strsep = $StrSep if (!defined $strsep);

  my ($angle, $str, $sign, $wholesec, $secfract, $min);

  if ($mode eq 'H') {
    $angle = $turn * 24;
  } else {
    $angle = $turn * 360;
  }

  if ($angle < 0.0) {
    $sign = -1;
    $angle = -$angle;
  } else {
    $sign = 1;
  }

  my $wholeangle = (int $angle);

  $angle -= $wholeangle;
  $angle *= 3600;

  # Get second fraction
  $wholesec = int $angle;
  $secfract = $angle - $wholesec;

  $wholesec %= 60;
  $min = ($angle-$wholesec - $secfract)/60.0;
  $secfract = int ($secfract * 10**$sig + 0.5); # Add 0.5 to ensure rounding

  # Check we have not rounded too far
  if ($secfract >= 10**$sig) {
    $secfract -= 10**$sig;
    $wholesec++;
    if ($wholesec >= 60.0) {
      $wholesec -= 60;
      $min++;
      if ($min >= 60.0) {
	$min -= 60;
	$wholeangle++;
      }
    }
  }

  my $angleform;
  if ($StrZero > 0) {
    $angleform = "%0$StrZero";
  } else {
    $angleform = '%';
  }

  my ($sep1, $sep2, $sep3);
  if ($strsep eq 'HMS') {
    if ($mode eq 'H') {
      $sep1 = 'H';
    } else {
      $sep1 = 'D';
    }
    $sep2 = 'M';
    $sep3 = 'S';
  } elsif ($strsep eq 'hms') {
    if ($mode eq 'H') {
      $sep1 = 'h';
    } else {
      $sep1 = 'd';
    }
    $sep2 = 'm';
    $sep3 = 's';
  } elsif ($strsep eq 'deg') { # What if $mode eq 'H'??
    $sep1 = 'd';
    $sep2 = "'";
    $sep3 = '"';
  } else {
    $sep1 = $sep2 = $strsep;
    $sep3 = '';
  }

  if ($sig > 0) {
    $str = sprintf("${angleform}d$sep1%02d".
		   "$sep2%02d.%0${sig}d$sep3", 
		   $wholeangle, $min, $wholesec, $secfract);
  } else {
    $str = sprintf("${angleform}d$sep1%02d".
		   "$sep2%02d$sep3", 
		   $wholeangle, $min, $wholesec);
  }

  if ($sign == -1) {
    $str = '-'.$str;
  }
  return $str;
}








sub mjd2cal($) {

  my $mjd = shift;

  my $ut = fmod($mjd,1.0);
  if ($ut<0.0) {
    $ut += 1.0;
    $mjd -= 1;
  }

  use integer;  # Calculations require integer division and modulation
  # Get the integral Julian Day number
  my $jd = nint($mjd + 2400001);

  # Do some rather cryptic calculations

  my $temp1 = 4*($jd+((6*(((4*$jd-17918)/146097)))/4+1)/2-37);
  my $temp2 = 10*((($temp1-237)%1461)/4)+5;

  my $year = $temp1/1461-4712;
  my $month =(($temp2/306+2)%12)+1;
  my $day = ($temp2%306)/10+1;

  return($day, $month, $year, $ut);
}






sub mpc_format_miti {
    # Format a MITI hash into MPC.
    my %stuff = @_;
    my ($day, $month, $year, $ut) = mjd2cal($stuff{EPOCH_MJD});   # cvt to cal date
    my $date = sprintf "%4d %02d %09.6f", $year, $month, $day + $ut;   # format YYYY MM DD.DDDDDD

    #my $saveZero = $Astro::Time::StrZero;   # save package global val
    my $saveZero = $StrZero;

    $StrZero = 2;  # 2 leading digits for stringified RA/DEC

    my $ra = deg2str($stuff{RA_DEG}, 'H', 3, ' ');
    print "Converted ra of $stuff{RA_DEG} degrees to $ra in HMS format\n";
    my $dec = ($stuff{DEC_DEG} > 0 ? '+' : '') . deg2str($stuff{DEC_DEG}, 'D', 2, ' ');

    my $line = sprintf
        $MPC_FORMAT_SPEC,
        "     ",                # number
        $stuff{ID},           # designation
        "",                     # discovery
        "",                     # code
        "C",                    # observation code
        $date,
        $ra,
        $dec,
        $stuff{MAG},
        "V",                    # MITI handles V only
        $stuff{OBSCODE},
        "";

    $StrZero = $saveZero;  # restore package global
    return $line;
}




sub mpc_format {
    # Format specified single observation table into single MPC-format line.
    # No tests are done to verify the data conforms to MPC format; namely
    # that the fields fit into allocated columns in the format.  If a
    # data item is larger than the columns allowed by the format, the data
    # item is truncated.
    my $obs = shift;
    print "obs = $obs \n";
    print " obs->RA = $obs->{RA} \n";
    exit(1);
    my $line = sprintf
        $MPC_FORMAT_SPEC,
        $obs->{NUMBER} ? sprintf("%05s", $obs->{NUMBER}) : "     ",
        $obs->{DESIGNATION},
        $obs->{DISCOVERY},
        $obs->{CODE},
        $obs->{OBSCODE},
        $obs->{DATE},
        $obs->{RA},
        $obs->{DEC},
        $obs->{MAG},
        $obs->{BAND},
        $obs->{OBSERVATORY};
    return $line;
}





sub mpc_write {
    # Write out list of observations to specified filename.
    # If filename is undef or '-', write to STDOUT.
    my ($stuff, $filename) = @_; 
    #my($stuff) = @_[0];
    #my($filename) = @_[1];
    print "stuff = $stuff, filename = $filename \n";
    #exit(1);
    my $str;
    my $line;

    if (not $filename or $filename eq '-') {
        foreach $line (@$stuff) {
            $str = mpc_format($line);    # format into single line
            print $str, "\n";
        }
    }
    else {
        open OUTFILE, ">$filename" or die "can't open $filename for writing";
        foreach $line (@$stuff) {
            $str = mpc_format($line);    # format into single line
            print OUTFILE $str, "\n";
        }
        close OUTFILE;
    }
}


use strict;

my($filename) = $ARGV[0];
my($outfile) = $ARGV[1];

my @dets;
my $count = 0;

open my $fh, $filename or die "can't open MITI file $filename";

open my $ofh, '>', $outfile or die "can't open outfile $outfile";

my $line;
my $miti;
my %mpc;
while (defined($line = <$fh>)) {

    if ($line !~ /^#/) {    # ignore comments
	$count = $count + 1;
	my(%miti) = miti_parse($line);
	
	#print "miti is $miti\n";
	#print "got det: MJD $miti->{EPOCH_MJD} RA $miti->{RA_DEG} DEC $miti->{DEC_DEG} \n";
	#while ( my ($key, $value) = each(%miti) ) {
	#    print "$key => $value\n";
	#}
	my($mpc) = mpc_format_miti(%miti);
	print $ofh $mpc . " \n";
	#exit(1);
	push @dets, \%mpc;
    }
}
print "count = $count\n";
close($ofh);

#print " dets has $#dets elements\n";
#exit(1);
#mpc_write(\@dets, $outfile)
