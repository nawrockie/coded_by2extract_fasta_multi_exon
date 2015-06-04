#!/opt/perl-5.16.2/bin/perl
#
# EPN, Fri Feb 27 14:52:13 2015
# 
# coded_by2extract_fasta_multi_exon.pl: expected input is output from
# edirect of the format: <protein-accession>
# <nucleotide-accession>:<start>-<stop> This script outputs much the
# same information in extract_fasta's expected input format.
# 
# Example input lines:
# AJJ35879.1	complement(CP009997.1:2173412..2176090)
# EEH43584.1	complement(join(KN275973.1:226623..226774, KN275973.1:226854..229725))
# AIJ19510.1	KJ818426.1:<1..>996
# YP_002836394.1	NC_012622.1:154158..154577
#
# (blank lines are skipped)
#
# Example output lines:
# CP009997.1 1 2173412 2176090 -
# KN275973.1 2 226623 226774 226854 229725 -
# KJ818426.1 1 1 996 +
# NC_012622.1 1 154158 154577 +

use strict;
use warnings;

my $line;          # a line of input, that we'll modify
my $orig_line;     # saved original line of input
my $strand;        # '+' or '-'
my $expected_accn; # accession of first exon, expected for all subsequent exons

while($line = <>) { 
  if($line =~ m/\w/) { 
    $orig_line = $line; # keep this saved so we can output it in case of an ERROR
    chomp $line;
    
    # Determine if 'complement' exists (which gives us strand information) and get 
    # rid of it if it does.
    # examples: 
    # AJJ35879.1	complement(CP009997.1:2173412..2176090)
    # EEH43584.1	complement(join(KN275973.1:226623..226774, KN275973.1:226854..229725))
    
    if($line =~ s/complement\(//) { # negative strand
      $strand = "-";
      if($line =~ s/\)$//) { ; } # remove final character, it should be a ')' matching the '(' in 'complement(' we just removed
      else { die "ERROR unable to parse (failure to find matching ) for 'complement ('; line: $orig_line"; }
    }
    else { # positive strand
      $strand = "+";
    }
    
    # Next deal with lines with 'join' indicating more than one subsequence (i.e. multiple exons)
    # example: 
    # EEH43584.1	join(KN275973.1:226623..226774, KN275973.1:226854..229725)
    
    if($line =~ s/^\S+\s+join\(//) { # get rid of everything up to 'join('
      # now: KN275973.1:226623..226774, KN275973.1:226854..229725))
      chomp $line;
      if($line =~ s/\)$//) { # remove trailing '))'
        # now: KN275973.1:226623..226774, KN275973.1:226854..229725
        # process each ', ' seperate element at a time
        my @elA = split(", ", $line);
        my $nel = scalar(@elA);
        my $start_stop = ""; # we'll append each exon to this 
        my $nexon = 0;
        my ($accn, $start, $stop);
        for(my $i = 0; $i < $nel; $i++) { 
          ($accn, $start, $stop) = breakdown($elA[$i], $orig_line);
          if($i == 0) { 
            $expected_accn = $accn; 
          }
          elsif($accn ne $expected_accn) {
            die "ERROR not all accessions in 'join' are identical in line: $orig_line";
          }
          $start_stop .= $start . " " . $stop . " ";
          $nexon++;
        }
        print "$accn $nexon $start_stop" . $strand . "\n"; 
      } # end of 'if($line =~ s/\)$//)' 
      else { 
        die "ERROR unable to parse (failure to find matching ) for 'join ('; line: $orig_line";
      }
    } # end of 'if($line =~ s/^\S+\s+\(join\(//) {'
    elsif($line =~ /^\S+\s+(.+)\:(\<?\d+)\.\.(\>?\d+)$/) { # no 'join' just one subsequence, easy case
      #AIV71043.1	CP009235.1:2698208..2701135
      $line =~ s/^\S+\s+//; 
      #now: CP009235.1:2698208..2701135
      my ($accn, $start, $stop) = breakdown($line, $orig_line);
      print "$accn 1 $start $stop $strand\n"; 
    }
    else { 
      die "ERROR unable to parse line: $orig_line";
    }
  }
}  

sub breakdown { 
# breakdown a token like this:
#
# CP009235.1:2698208..2701135
# JGZA01000019.1:<1..1949 # flush with beginning of sequence
# GQ358600.1:<1..>453     # flush with end of sequence
# 
# into just <accn> <start> <stop>
# and determine if the subsequence ends (start and/or stop) are
# flush with the beginning or end of the source sequence.
# as indicated by a '<' before the start coordinate and 
# a '>' before the end coordinate.
#
# Return: 
# $accn:       the accession
# $start:      the start coordinate
# $stop:       the stop coordinate
# $flushbegin: '1' if start is first nt of source sequence, else '0'
# $flushend:   '1' if start is first nt of source sequence, else '0'
#
# Return values for examples above:
# 
# CP009235.1 2698208 2701135 0 0 
# JGZA01000019.1 1 1949 1 0 
# GQ358600.1:1 453 1 1
# 
  if(scalar(@_) != 2) { die "ERROR entered breakdown() with wrong number of args"; }
  my ($token, $line) = (@_);
  my ($flushbegin, $flushend);
  
  if($token =~ /^(.+)\:(\<?\d+)\.\.(\>?\d+)/) { 
    my ($accn, $start, $stop) = ($1, $2, $3);
    if($start =~ s/^\<//) { $flushbegin = 1; }
    else                  { $flushbegin = 0; }
    if($stop  =~ s/^\>//) { $flushend = 1; }
    else                  { $flushend = 0; }
    if($start > $stop) { die "ERROR unexpectedly start > stop in $token on line: $orig_line"; }
    return ($accn, $start, $stop, $flushbegin, $flushend);
  }
  else { 
    die "ERROR unable to breakdown $token, in line: $line"; 
  }
}
  
