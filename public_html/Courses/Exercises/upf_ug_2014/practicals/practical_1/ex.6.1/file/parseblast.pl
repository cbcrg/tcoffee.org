#!/usr/bin/perl -w
# 
# ##################################################################
# #                          parseblast                            #
# ##################################################################
# 
#             Extracting HSPs from blast programs output.
# 
#       Copyright (C) 2000-2003, Josep Francesc ABRIL FERRANDO  
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# ##################################################################
# 
# $Id: parseblast.nw,v 1.10 2004/05/15 00:02:12 jabril Exp $
#
use strict;

#
#
use Getopt::Long;
Getopt::Long::Configure("bundling");

#
use Benchmark;
my @exectime = ();
&init_timer( \@exectime );

##############################################################################
##                              GLOBAL VARS                                 ##
##############################################################################
#
my $PROGRAM = "parseblast";
my $VERSION = "v2.1";
my $PLVER = sprintf( "v%vd", $^V );
my $DATE  = localtime;
my $CDATE = &get_date();
my ( $USER, $HOST );
if ( defined( $ENV{USER} ) ) {
    $USER = $ENV{USER};
}
else {
    chomp( $USER = system("whoami") );
}
if ( defined( $ENV{HOSTNAME} ) ) {
    $HOST = $ENV{HOSTNAME};
}
else {
    chomp( $HOST = system("hostname") );
}

#
# FLAGS
my (
    $hsp_flg,       $gff_flg,      $fullgff_flg,  $aplot_flg,
    $nogff_flg,     $subject_flg,  $sequence_flg, $full_scores,
    $comment_flg,   $nocmmnt_flg,  $split_flg,    $help_flg,
    $ver_flg,       $err_flg,      $hitbl_flg,    $expanded_flg,
    $short_tgt_flg, $pairwise_flg, $msf_flg,      $aln_flg,
    $bit_flg,       $ids_flg,      $noframe_flg
);
my ( $prt, $main, $seqflg, $hsp, $fragment, $param, $aln_split, $prt_pos_flg ) =
  ( 0, 0, 0, 0, 0, 0, 0, 0 );

#
# Parser Vars
my ( $ifh, $sin_flg, $file );
my @ifiles = ();
my ( $fs,        $GFFfmt,      $GFFtgs,  $GFFtgl,  $GFFtg );
my ( $prt_funct, $prog_params, $program, $version, $seqname );
my ( $scQ, $scS );    # reSCale_Query/Subject : 1 to re-scalate lengths.
my ( $query_name, $db_name, $score, $descr, $ori, $end, $seq, $txt, $tt, $CT,
    $PT, $HT );
my (
    @seqlist, %prgseq, %dbase,     %query,   %cnt,
    %desc,    %sco,    %hsp_start, %hsp_end, %hsp_seq
);
my ( $qm, $sm, $x, $y, $ml, $a, $b, $aq, $as, $sql, $lq, $ls, $sq );
my $foe            = 0;
my $index          = 0;
my $chars_per_line =
  50;    # chars of sequences to show per line (in alignment modes)
my ( $tlen, $clen ) = ( 54, 50 );

#
# Formatting output as plain, HSPs, GFF, APLOT, ALN.
my (
    $n,  $nm, $wnm, $tq, $ts,    # couNter, NaMe, TagQuery, TagSubject
    $sc, $bt, $ex,  $pv, $id,    # SCore, BiTscore, EvalueX, IDentityscore
    $scores,                     # showing all scores
    $frq, $frs,                  # QueryFRame, SubjectFRame
    $stq, $sts,                  # QuerySTrand, SubjectSTrand
    $gsc, $prg,                  # GroupSCore, PRoGram
    $hsq, $hss, $heq,
    $hes,    # HspStartQuery, HspStartSubject, HspEndQuery, HspEndSubject
    $lnq, $lns, $lnmin, $lnmax
    , # LeNgthQuery, LeNgthSubject, LeNgthMINqueyxsubject, LeNgthMAXqueyxsubject
    $lnmx, $gpq, $gps, $gpt    # LeNgthMaXhspseq, GaPQuery, GaPSubject, GaPTotal
);

##############################################################################
##                              MAIN LOOP                                   ##
##############################################################################
#
&init_main_vars();
&get_commandline_opts();    # if help is ON no start header must be shown...
&program_started("$PROGRAM $VERSION");
foreach $file (@ifiles) {
    ( $ifh, $sin_flg ) = &getIFH($file);
    &parseblast( $ifh, $file );
    close($ifh) unless $sin_flg == 2;
}
&program_finished("$PROGRAM $VERSION");
exit(0);

##############################################################################
##                       Getting Comand-line Options                        ##
##############################################################################
#
sub get_commandline_opts() {
    GetOptions(
        "G|gff"            => \$gff_flg,
        "F|fullgff"        => \$fullgff_flg,
        "A|aplot"          => \$aplot_flg,
        "S|subject"        => \$subject_flg,
        "Q|sequence"       => \$sequence_flg,
        "H|hit-table"      => \$hitbl_flg,
        "X|extended"       => \$expanded_flg,
        "P|pairwise"       => \$pairwise_flg,
        "M|msf"            => \$msf_flg,
        "N|aln"            => \$aln_flg,
        "W|show-coords"    => \$prt_pos_flg,
        "b|bit-score"      => \$bit_flg,
        "i|identity-score" => \$ids_flg,
        "s|full-scores"    => \$full_scores,
        "u|no-frame"       => \$noframe_flg,
        "t|compact-tags"   => \$short_tgt_flg,
        "c|comments"       => \$comment_flg,
        "n|no-comments"    => \$nocmmnt_flg,
        "v|verbose"        => \$err_flg,
        "version"          => \$ver_flg,
        "h|help|\?"        => \$help_flg
    );

    #
    ($ver_flg)  && &prt_version;
    ($help_flg) && &prt_help;

    #
    FLGS: {    # first choose disables any other command-line option.
        $aln_flg
          && ( $expanded_flg = $msf_flg = $pairwise_flg = $hsp_flg =
            $comment_flg = 0, $nogff_flg = 1, last FLGS );
        $msf_flg
          && ( $expanded_flg = $aln_flg = $pairwise_flg = $hsp_flg =
            $comment_flg = 0, $nogff_flg = 1, last FLGS );
        $pairwise_flg
          && ( $expanded_flg = $aln_flg = $msf_flg = $hsp_flg = $comment_flg =
            0, $nogff_flg = 1, last FLGS );
        $gff_flg
          && ( $expanded_flg = $aln_flg = $msf_flg = $pairwise_flg =
            $fullgff_flg = $aplot_flg = $hsp_flg = 0, $nogff_flg = 0,
            last FLGS );
        $fullgff_flg
          && ( $expanded_flg = $aln_flg = $msf_flg = $pairwise_flg = $gff_flg =
            $aplot_flg = $hsp_flg = 0, $nogff_flg = 0, last FLGS );
        $aplot_flg
          && ( $expanded_flg = $aln_flg = $msf_flg = $pairwise_flg = $gff_flg =
            $fullgff_flg = $hsp_flg = 0, $nogff_flg = 0, last FLGS );
        $expanded_flg
          && ( $aln_flg = $msf_flg = $pairwise_flg = $hsp_flg = 0,
            $nogff_flg = 1, last FLGS );

        # $hitbl_flg => $hsp_flg
        $hsp_flg   = 1;
        $nogff_flg = 1;
    };    # FLGS
          #
    FLGS: {    # first choose disables any other command-line option.
        $bit_flg && ( $ids_flg = 0, last FLGS );
        $ids_flg && ( $bit_flg = 0 );
    };

    #    
    $fs = $fullgff_flg ? "\t" : " ";
    $GFFfmt = ( ("\%s${fs}") x 8 ) . "\%s" . "\n";
    $GFFtgs = " \%s \%s \%s \%s;${fs}";
    $GFFtgl =
      ";${fs}Start \%s;${fs}End \%s;${fs}Strand \%s;${fs}Frame \%s;${fs}";
    $GFFtg = $short_tgt_flg ? $GFFtgs : $GFFtgl;

    #
    $nocmmnt_flg && ( $comment_flg = 0 );

    # 
    @ifiles = ( scalar(@ARGV) > 0 ? @ARGV : "-" );
}    # get_commandline_opts

##############################################################################
##                          Global Subroutines                              ##
##############################################################################
#
sub program_started() {
    my $prog = shift;
    &header(
        "RUNNING $prog", '', "Host: $HOST", "User: $USER",
        "Perl: $PLVER",  '', "Date: $DATE"
      )
      if $err_flg;
}    # program_started

sub program_finished() {
    my $prog = shift;
    my $txt = &timing( \@exectime, 1 );
    $txt =~ s{^\s*(\d+)\s+wallclock\s+secs\s+(\()}{\n$2}o && do {
        $txt = sprintf(
            "Total Execution Time: %02d:%02d:%02d (%d secs)%s",
            sub { ( $_[2] - 1, $_[1], $_[0] ) }
            ->( localtime($1) ), $1, $txt
        );
    };
    &header( "$prog HAS FINISHED", '', ( split /\n/o, $txt ) ) if $err_flg;
}    # program_finished
     #

sub header() {
    my ( @lns, $str, $t, $comment );
    @lns     = @_;
    $str     = ( "#" x ( $tlen + 4 ) ) . "\n";
    $t       = ( "#" x 2 );
    $comment = $str;
    foreach my $ln (@lns) {
        $comment .= "$t" . ( &fill_mid( $ln, $tlen, " " ) ) . "$t\n";
    }
    print STDERR $comment . $str;
}    # header
     #
sub fill_right { $_[0] . ( $_[2] x ( $_[1] - length( $_[0] ) ) ) }
sub fill_left { ( $_[2] x ( $_[1] - length( $_[0] ) ) ) . $_[0] }

sub fill_mid() {
    my $l = length( $_[0] );
    my $k = int( ( $_[1] - $l ) / 2 );
    return ( $_[2] x $k ) . $_[0] . ( $_[2] x ( $_[1] - ( $l + $k ) ) );
}    # fill_mid
     #

sub prt_progress {
    $err_flg && do {
        print STDERR ".";
        ( ( $_[0] % $clen ) == 0 )
          && print STDERR "[" . &fill_left( $_[0], 6, "0" ) . "]\n";
    };
}    # prt_progress

sub prt_foeprg {
    $err_flg && do {
        ( ( $_[0] % $clen ) != 0 )
          && print STDERR "[" . &fill_left( $_[0], 6, "0" ) . "]\n";
    };
}    # prt_foeprg
     #

sub max {
    my ($z) = shift @_;
    my $l;
    foreach $l (@_) { $z = $l if $l > $z; }
    return $z;
}    # max

#
# Timing.
sub init_timer() {
    my $refary = shift;
    @{$refary} = ( new Benchmark );
    return;
}    # init_timer

sub timing() {
    my ( $refary, $tmp ) = @_;
    my $flg = defined($tmp) || 0;
    push @{$refary}, ( new Benchmark );
    my $mx = $#{$refary};

    # partial time
    $flg || do {
        return timestr( timediff( $refary->[$mx], $refary->[ ( $mx - 1 ) ] ) );
    };

    # total time
    return timestr( timediff( $refary->[$mx], $refary->[0] ) );
}    # timing

sub get_date() {
    return sprintf(
        "%04d%02d%02d%02d%02d%02d", sub {
            $_[5] + 1900, $_[4] + 1, $_[3], $_[2], $_[1], $_[0];
        }
        ->(localtime)
    );
}    # get_date

#
# Print help to STDERR.
sub prt_version() {
    print STDERR "#\n# $PROGRAM $VERSION \n#\n";
    exit(2);
}    # prt_version

sub prt_help() {
    open( HELP, "| more" );
    print HELP <<EndOfHelp;
PROGRAM:
        $PROGRAM $VERSION

USAGE:  parseblast.pl [options] <results.from.blast>

DESCRIPTION:

  Filtering High-scoring Segment Pairs (HSPs) from WU/NCBI BLAST.
  Different output options are available, the most important here are
  those allowing to write HSPs in GFF format (GFFv1, GFFv2 or APLOT).
  Sequences can be included in the GFF records as a comment field.
  Furthermore, this script can output also the alignments
  for each HSP in ALN, MSF or tabular formats. 

  NOTE.- If first line from blast program output (the one containing
      which flavour has been run, say here BLASTN, BLASTP, BLASTX,
      TBLASTN or TBLASTX), is missing, the program assumes that it
      contains BLASTN HSP records. So that, ensure that you feed
      the $PROGRAM script with a well formatted BLAST file.

         Sometimes there are no spaces between the HSP coords and
      its sequence, as it sometimes happens in Web-Blast or
      Paracel-Blast outputs. Now those records are processed ok
      and that HSP is retrieved as well as "standard" ones.

  WARNING.- Frame fields from GFF records generated with $PROGRAM
      contain BLAST frame (".","1","2","3") instead of the GFF standard
      values (".","0","1","2"). As the frame for reverse strand must be
      recalculated from the original sequence length, we suggest
      users to post-process the GFF output from this script with
      a suitable filter that fix the frames (in case that the program
      that is going to use the GFF records will not work with 
      the original BLAST frames). We provide the command-line option
      "--no-frame" to set frames to "." (meaning that there is no frame).

COMMAND-LINE OPTIONS:

    $PROGRAM prints output in "HSP" format by default (see below).
  It takes input from <STDIN> or single/multiple files, and writes
  its output to <STDOUT>, so user can redirect to a file but
  he also could use the program as a filter within a pipe. 
    "-N", "-M", "-P", "-G", "-F", "-A" and "-X" options (also the long
  name versions for each one) are mutually exclusive, and their
  precedence order is shown above.

  GFF OPTIONS:
    -G, --gff            : prints output in GFFv1 format.
    -F, --fullgff        : prints output in GFFv2 "alignment" format ("target").
    -A, --aplot          : prints output in pseudo-GFF APLOT "alignment" format.
    -S, --subject        : projecting GFF output by SUBJECT (default by QUERY).
    -Q, --sequence       : append query and subject sequences to GFF record.
    -b, --bit-score      : set <score> field to Bits (default Alignment Score).
    -i, --identity-score : set <score> field to Identities (default Alignment).
    -s, --full-scores    : include all scores for each HSP in each GFF record.
    -u, --no-frame       : set all frames to "." (GFF for not available frames).
    -t, --compact-tags   : target coords+strand+frame in short form (NO GFFv2!).

  ALIGNMENT OPTIONS:
    -P, --pairwise       : prints pairwise alignment for each HSP in TBL format.
    -M, --msf            : prints pairwise alignment for each HSP in MSF format.
    -N, --aln            : prints pairwise alignment for each HSP in ALN format.
    -W, --show-coords    : adds start/end positions to alignment output.

  GENERAL OPTIONS:
    -H, --hit-table      : "HITBL" output, similar to blast hit table output
                           but simpler than default "HSP" format.
    -X, --expanded       : "EXPANDED" output (generating multiline records).
    -c, --comments       : include parameters from blast program as comments.
    -n, --no-comments    : do not print "#" lines (raw output without comments).
    -v, --verbose        : warnings sent to <STDERR>.
        --version        : prints program version and exits.
    -h, --help           : shows this help and exits.

OUTPUT FORMATS:

    "S_" stands for "Subject_Sequence" and "Q_" for "Query_Sequence". <Program>
  name is taken from input blast file. <Strands> are calculated from <start> and
  <end> positions on original blast file. <Frame> is obtained from the blast 
  file if is present else is set to ".". <SCORE> is set to Alignment Score by 
  default, you can change it with "-b" and "-i".
    If "-S" or "--subject" options are given, then QUERY fields are referred to
  SUBJECT and SUBJECT fields are relative to QUERY (this only available for GFF
  output records).
    Dots ("...") mean that record description continues in the following line,
  but such record is printed as a single line record by $PROGRAM.

[HSP]  <- (This is the DEFAULT OUTPUT FORMAT)
 <Program> <DataBase> : ...
   ... <IdentityMatches> <Min_Length> <IdentityScore> ...
   ... <AlignmentScore> <BitScore> <E_Value> <P_Sum> : ...
   ... <Q_Name> <Q_Start> <Q_End> <Q_Strand> <Q_Frame> : ...
   ... <S_Name> <S_Start> <S_End> <S_Strand> <S_Frame> : <S_FullDescription>

[GFF]
 <Q_Name> <Program> hsp <Q_Start> <Q_End> <SCORE> <Q_Strand> <Q_Frame> <S_Name>

[FULL GFF]  <- (GFF showing alignment data)
 <Q_Name> <Program> hsp <Q_Start> <Q_End> <SCORE> <Q_Strand> <Q_Frame> ...
   ... Target "<S_Name>" <S_Start> <S_End> ...
   ... E_value <E_Value> Strand <S_Strand> Frame <S_Frame>

[APLOT]  <- (GFF format enhanced for APLOT program)
 <Q_Name>:<S_Name> <Program> hsp <Q_Start>:<S_Start> <Q_End>:<S_End> <SCORE> ...
   ... <Q_Strand>:<S_Strand> <Q_Frame>:<S_Frame> <BitScore>:<HSP_Number> ...
   ... \# E_value <E_Value>

[HITBL]
 <Q_Name> <Q_Start> <Q_End> <Q_Strand> <Q_Frame> ...
   ... <S_Name> <S_Start> <S_End> <S_Strand> <S_Frame> ...
   ... <IdentityScore> <AlignmentScore> <BitScore> <E_Value> <P_Sum>

[EXPANDED]
 MATCH(<HSP_Number>): <Q_Name> x <S_Name>
 SCORE(<HSP_Number>): <AlignmentScore>
 BITSC(<HSP_Number>): <BitScore>
 EXPEC(<HSP_Number>): <E_Value> Psum(<P_Sum>)
 IDENT(<HSP_Number>): <IdentityMatches>/<Min_Length> : <IdentityScore> \%
 T_GAP(<HSP_Number>): <TotalGaps(BothSeqs)>
 FRAME(<HSP_Number>): <Q_Frame>/<S_Frame>
 STRND(<HSP_Number>): <Q_Strand>/<S_Strand>
 MXLEN(<HSP_Number>): <Max_Length>
 QUERY(<HSP_Number>): length <Q_Length> : gaps <Q_TotalGaps> : ...
   ... <Q_Start> <Q_End> : <Q_Strand> : <Q_Frame> : <Q_FullSequence>
 SBJCT(<HSP_Number>): length <S_Length> : gaps <S_TotalGaps> : ...
   ... <S_Start> <S_End> : <S_Strand> : <S_Frame> : <S_FullSequence>

BUGS:    Report any problem to: abril\@imim.es

AUTHOR:  $PROGRAM is under GNU-GPL (C) 2000-2003, Josep F. Abril

EndOfHelp
    close(HELP);
    exit(2);
}    # prt_help
 # print STDOUT "HSP:$hsp_flg GFF:$gff_flg COM:$comment_flg SPL:$split_flg HLP:$help_flg\n";

#
# Retrieving file handles
sub getIFH() {
    my $file = shift;
    my ( $fstr, $ferr, $sflg );
    local *D;
    defined($file) || do {
        &header(" WARNING!!! FILE NAME NOT DEFINED !!!\n");
        exit(1);
    };
    ( $file ne '-' ) || do {
        &header("## Input Stream set to STDIN...") if $err_flg;
        return ( *STDIN, 2 );
    };
    if ( $file =~ /\|\s*$/o ) {
        $fstr = $file;
        $ferr = 'NOREADPIPE';
        $sflg = 1;
    }
    elsif ( !-e $file ) {
        &header(" WARNING!!! INPUT File does not exist: \"$file\"\n");
        exit(1);
    }
    elsif ( $file =~ /\.gz\s*$/o ) {
        $fstr = "zcat $file |";
        $ferr = 'NOREADPIPE';
        $sflg = 1;
    }
    else {
        $fstr = "< $file";
        $ferr = 'NOREADFILE';
        $sflg = 0;
    }
    open( D, $fstr ) || do {
        &header(" WARNING!!! Cannot open input stream \"$fstr\"");
        exit(1);
    };
    $err_flg && print STDERR "##\n##\n##Input Stream set to: \"$fstr\"\n##\n";
    return ( *D, $sflg );
}    # getIFH

#
# Get new lines while empty line is not found and append to last line.
sub get_lines {
    my ( $spc, $ifh ) = @_;
    my $tmp;

    # local($tmp);
    while (<$ifh>) {
        last if /^\s*$/;

        # print STDERR "$_";
        chomp;
        $_ =~ s/^\s*/$spc/og;
        $tmp .= $_;
    }
    return $tmp;
}    # get_lines

#
# Getting scores from scoring vector extracted from HSP record.
sub get_scores {
    my $t = $_[0];
    my ( $sc, $bt, $ex, $pv, $id );
    my ( $qfr, $sfr ) = ( '.', '.' );
    ( ( $t =~ /Score[^\s]*\s+=\s+\b(\d+)\b\s+\([^,]*,/ )
      || ( $t =~ /Score\s+=\s+[^,]*\s+\((\d+)\)[^,]*,/o ) ) && ( $sc = $1 );
    ( $t =~ /Score[^\s]*\s+=.*[\s\(]([+-]?(\d+\.?\d*|\.\d+))\b \bbits[^,]*,/o )
      && ( $bt = $1 );
    ( $t =~
/Expect[^\s]*\s+=\s+([+-]?([Ee][+-]?\d+|(\d+\.?\d*|\.\d+)([Ee][+-]?\d+)?))\s*,/o
    ) && ( $ex = $1 );
    ( $ex =~ /^[Ee]/o )    && ( $ex = "1" . $ex );
    ( $ex =~ /^(.*\.)$/o ) && ( $ex = $1 . "0" );
    ( $t  =~
/Sum[^\s]*\s+\bP[^\s]*\s+=\s+([+-]?([Ee][+-]?\d+|(\d+\.?\d*|\.\d+)([Ee][+-]?\d+)?))\s*,/o
    ) ? ( $pv = $1 ) : ( $pv = "." );    # $pv not defined then $pv="."
    ( $pv =~ /^(.*\.)$/o ) && ( $pv = $1 . "0" );
    ( $t =~ /Identities[^\s]*\s+=\s+(\d+)\/(\d+)\s+/o ) && ( $id = $1 );
    $noframe_flg || do {

        # BLASTX (translated nucleotides vs protein)
        ( $scQ && !$scS ) && do {
            ( $t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\b/o ) && ( $qfr = $2 );
        };

        # TBLASTN (protein vs translated nucleotides)
        ( !$scQ && $scS ) && do {
            ( $t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\b/o ) && ( $sfr = $2 );
        };

        # TBLASTX (translated nucleotides vs translated nucleotides)
        ( $scQ && $scS ) && do {
            ( $t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\s*\/\s*(\+|\-)(\d)\b/o )
              && ( $qfr = $2, $sfr = $4 );
        };
    };    # $noframe_flg
    return ( $sc, $bt, $ex, $pv, $id, $qfr, $sfr );
}    # get_scores

#
# Calculating strand from coords
sub chk_strand {
    my ( $first, $last ) = @_;
    my $st = "+";
    ( $first >= $last ) && ( ( $first, $last ) = ( $last, $first ), $st = "-" );
    return ( $first, $last, $st );
}    # chk_strand

#
# Formatting output as plain, HSPs, GFF, APLOT, ALN.
sub prt_hsp {
    print STDOUT <<"EndOfHSPs";
$prg $dbase{$nm} : $id $lnmin $gsc $sc $bt $ex $pv : $query{$nm} $hsq $heq $stq $frq : $wnm $hss $hes $sts $frs : $desc{$nm}
EndOfHSPs
}    # prt_hsp
     #

sub prt_hittbl {
    print STDOUT <<"EndOfHSPs";
$query{$nm} $hsq $heq $stq $frq $wnm $hss $hes $sts $frs $gsc $sc $bt $ex $pv
EndOfHSPs

    # $id $lnmin $lnmx $lnq $lns $gpt $gpq $gps
}    # prt_hittbl
     #

sub prt_ext {
    print STDOUT <<"EndOfPlain";
MATCH($n): $query{$nm} x $wnm
SCORE($n): $sc\nBITSC($n): $bt\nEXPEC($n): $ex Psum($pv)
IDENT($n): $id/$lnmin : $gsc \%
T_GAP($n): $gpt\nFRAME($n): $frq/$frs\nSTRND($n): $stq/$sts\nMXLEN($n): $lnmx
QUERY($n): length $lnq : gaps $gpq : $hsq $heq : $stq : $frq : $hsp_seq{$tq}
SBJCT($n): length $lns : gaps $gps : $hss $hes : $sts : $frs : $hsp_seq{$ts}
\#\#
EndOfPlain
}    # prt_ext
     #

sub prt_Q_gff {    # $hss $hes;\tStrand $sts;\tFrame $frs;\t
    printf STDOUT $GFFfmt, $query{$nm}, $prg, 'hsp', $hsq, $heq, $gsc, $stq,
      $frq,
      "$wnm" . ( sprintf( $GFFtg, $hss, $hes, $sts, $frs ) ) . "$scores$sq";
}

sub prt_S_gff {    # $hsq $heq;\tStrand $stq;\tFrame $frq;\t
    printf STDOUT $GFFfmt, $wnm, $prg, 'hsp', $hss, $hes, $gsc, $sts, $frs,
      "$query{$nm}"
      . ( sprintf( $GFFtg, $hsq, $heq, $stq, $frq ) )
      . "$scores$sq";
}

sub prt_Q_fullgff {    # $hss $hes;\tStrand $sts;\tFrame $frs;\t
    printf STDOUT $GFFfmt, $query{$nm}, $prg, 'hsp', $hsq, $heq, $gsc, $stq,
      $frq,
      "Target \"$wnm\""
      . ( sprintf( $GFFtg, $hss, $hes, $sts, $frs ) )
      . "$scores$sq";
}

sub prt_S_fullgff {    # $hsq $heq;\tStrand $stq;\tFrame $frq;\t
    printf STDOUT $GFFfmt, $wnm, $prg, 'hsp', $hss, $hes, $gsc, $sts, $frs,
      "Target \"$query{$nm}\""
      . ( sprintf( $GFFtg, $hsq, $heq, $stq, $frq ) )
      . "$scores$sq";
}

sub prt_Q_aplot {
    printf STDOUT $GFFfmt, "$query{$nm}:$wnm", $prg, 'hsp', "$hsq:$hss",
      "$heq:$hes", $gsc, "$stq:$sts", "$frq:$frs", "$bt:$n $scores$sq";
}

sub prt_S_aplot {
    printf STDOUT $GFFfmt, "$wnm:$query{$nm}", $prg, 'hsp', "$hss:$hsq",
      "$hes:$heq", $gsc, "$sts:$stq", "$frs:$frq", "$bt:$n $scores$sq";
}

sub prt_pairwise {
    $prt_pos_flg && do {
        $ml = &max( length($hsq), length($heq), length($hss), length($hes) );
        ( $x, $y ) =
          ( " "
          . &fill_left( $hsq, $ml, " " ) . " "
          . &fill_left( $heq, $ml, " " )
          . " $stq $frq ",
          " "
          . &fill_left( $hss, $ml, " " ) . " "
          . &fill_left( $hes, $ml, " " )
          . " $sts $frs " );
    };
    print "#\n" if $aln_split;
    $aln_split = 1;
    print STDOUT <<"EndOfALIGN";
$a$x$hsp_seq{$tq}
$b$y$hsp_seq{$ts}
EndOfALIGN
}    # prt_pairwise
     #

sub prt_cmn_aln {
    print STDOUT <<"EndOfALN" if !$nocmmnt_flg;
################################################################################
##
##  $query{$nm} x $wnm #$n
##
##    $prg $dbase{$nm}
##    Identity: $gsc  Score: $sc  Bits: $bt  E_value: $ex  P_value: $pv
##    DESCR: $desc{$nm}
##
EndOfALN
}    # prt_cmn_aln

sub prt_msf {
    my ( $i, $jq, $js, $vs, $ve, $chw );
    ( $aq, $as, $sql ) =
      ( $hsp_seq{$tq}, $hsp_seq{$ts}, length( $hsp_seq{$tq} ) );
    $aq =~ s/\.//og;
    $as =~ s/\.//og;
    ( $lq, $ls ) = ( length($aq), length($as) );
    &prt_cmn_aln;
    print STDOUT <<"EndOfMSF";
\n
$query{$nm}\_x_$wnm\_\#$n.msf  MSF: $sql  Type: P  May 4th, 2000  Check: 0  ..\n
Name: $a  Len: $lq  Check: 0  Weight: 1.00
Name: $b  Len: $ls  Check: 0  Weight: 1.00
\n//
EndOfMSF
    for ( $i = 0 ; $i <= ( $sql - 1 ) ; $i += $chars_per_line ) {
        ( $jq, $js ) =
          ( substr( $hsp_seq{$tq}, $i, 50 ), substr( $hsp_seq{$ts}, $i, 50 ) );
        print STDOUT "\n";
        $prt_pos_flg && do {
            $chw = length($jq);
            ( $vs, $ve ) = ( $i + 1, $i + $chw );
            $chw = $chw - length( $vs . $ve );
            print STDOUT " " x ( $ml + 1 ) . $vs
              . " " x ( $chw > 1 ? $chw : 1 ) . $ve . "\n";
        };
        print STDOUT <<"EndOfALIGN";
$a $jq
$b $js
EndOfALIGN
    }
    print STDOUT "\n\n";
}    # prt_msf

sub prt_aln {
    my ( $i, $jq, $js, $vs, $ve, $chw );
    ( $aq, $as, $sql ) =
      ( $hsp_seq{$tq}, $hsp_seq{$ts}, length( $hsp_seq{$tq} ) );
    $aq =~ s/\./-/og;
    $as =~ s/\./-/og;
    &prt_cmn_aln;
    for ( $i = 0 ; $i <= ( $sql - 1 ) ; $i += $chars_per_line ) {
        ( $jq, $js ) = ( substr( $aq, $i, 50 ), substr( $as, $i, 50 ) );
        print STDOUT "\n";
        $prt_pos_flg && do {
            $chw = length($jq);
            ( $vs, $ve ) = ( $i + 1, $i + $chw );
            $chw = $chw - length( $vs . $ve );
            print STDOUT " " x ( $ml + 1 ) . $vs
              . " " x ( $chw > 1 ? $chw : 1 ) . $ve . "\n";
        };
        print STDOUT <<"EndOfALIGN";
$a $jq
$b $js
EndOfALIGN
    }
    print STDOUT "\n\n";
}    # prt_aln
     #

sub setprtfunct() {
    CHOOSE: {
        if ($hsp_flg) {
            $prt_funct =
              $hitbl_flg ? sub { &prt_hittbl() } : sub { &prt_hsp() };
        }
        elsif ( !$nogff_flg ) {
            if ($subject_flg) {
                $gff_flg && ( $prt_funct = sub { &prt_S_gff() }, last CHOOSE );
                $fullgff_flg
                  && ( $prt_funct = sub { &prt_S_fullgff() }, last CHOOSE );
                $aplot_flg
                  && ( $prt_funct = sub { &prt_S_aplot() }, last CHOOSE );
            }
            else {
                $gff_flg && ( $prt_funct = sub { &prt_Q_gff() }, last CHOOSE );
                $fullgff_flg
                  && ( $prt_funct = sub { &prt_Q_fullgff() }, last CHOOSE );
                $aplot_flg
                  && ( $prt_funct = sub { &prt_Q_aplot() }, last CHOOSE );
            }
        }
        elsif ($expanded_flg) {
            $prt_funct = sub { &prt_ext() };
        }
        else {
            $pairwise_flg
              && ( $prt_funct = sub { &prt_pairwise() }, last CHOOSE );
            $msf_flg && ( $prt_funct = sub { &prt_msf() }, last CHOOSE );
            $aln_flg && ( $prt_funct = sub { &prt_aln() }, last CHOOSE );
        }
    };    # CHOOSE
}    # setprtfunct
     #

sub prt_out() {
    my $ht = 0;
    $CT++;
    $prt = 0;
    &header("WRITING OUTPUT TO STDOUT") if $err_flg;
    &setprtfunct();
    while (@seqlist) {
        $nm = shift (@seqlist);
        ( $wnm = $nm ) =~ s/\_\d+$//o;
        ( !$hsp_flg && $comment_flg ) && do {
            print STDOUT "#\n# $prgseq{$nm} :: DB $dbase{$nm} :: $cnt{$nm} "
              . "HSPs for $query{$nm}x$wnm \n# DESCR: $desc{$nm}\n#\n";
        };
        ( $cnt{$nm} > 0 ) && do {
            for ( $n = 1 ; $n <= $cnt{$nm} ; $n++ ) {
                &prt_progress( ++$ht );
                $tq = $nm . "query" . $n;
                $ts = $nm . "sbjct" . $n;
                ( $sc, $bt, $ex, $pv, $id, $frq, $frs ) =
                  &get_scores( $sco{ $nm . $n } );
                ( $hsq, $heq, $stq ) =
                  &chk_strand( $hsp_start{$tq}, $hsp_end{$tq} );
                ( $hss, $hes, $sts ) =
                  &chk_strand( $hsp_start{$ts}, $hsp_end{$ts} );
                $lnq = $heq - $hsq + 1;
                $scQ && ( $lnq = $lnq / 3 );
                $lns = $hes - $hss + 1;
                $scS && ( $lns = $lns / 3 );
                $lnmin = ( $lnq > $lns ) ? $lns : $lnq;
                $lnmax = ( $lnq < $lns ) ? $lns : $lnq;
                $lnmx = length( $hsp_seq{$tq} );
                {
                    my $hh = $hsp_seq{$tq};
                    $gpq = ( $hh =~ s/-/ /og ) || 0;
                    $hh = $hsp_seq{$ts};
                    $gps = ( $hh =~ s/-/ /og ) || 0;
                };
                $gpt = $gpq + $gps;
                {
                    ( $ids_flg || $expanded_flg || $hsp_flg )
                      &&    # score is Identities divided by minlength
                      ( ($gsc) =
                        eval( ( $id / $lnmin ) * 100 ) =~ /^(\d+(\.\d{0,3})?)/o,
                        last );
                    $bit_flg && ( $gsc = $bt, last );
                    $gsc = $sc;
                };
                ($prg) = $prgseq{$nm} =~ /^([^\s]+)\s/o;

                # GFF format
                PRINT: {
                    $hsp_flg
                      && do { $prt_funct->(); last PRINT; };    # default output
                                                                #
                    do {
                        $sq = $scores = "";
                        $full_scores && do {
                            $scores =
                              "E_value $ex;${fs}P_sum $pv;${fs}Aln_Score"
                              . " $sc;${fs}Bit_Score $bt;${fs}Idn_Score "
                              . sprintf( "%6.2f", ( $id / $lnmin ) * 100 )
                              . " ($id/$lnmin);${fs}Gaps "
                              . sprintf( "%4.2f", ( $gpt / $lnmx ) * 100 )
                              . " (Q:$gpq|S:$gps);"
                              . "${fs}Lengths Q:$lnq|S:$lns|T:$lnmx";
                            $aplot_flg && ( $scores = "# $scores" );
                        };    # $full_scores
                        $subject_flg && do {
                            $sq =
                              ";${fs}Subject_seq $hsp_seq{$ts} ;${fs}"
                              . "Query_seq $hsp_seq{$tq}"
                              if $sequence_flg;
                            $prt_funct->();
                            last PRINT;
                        };    # $subject_flg
                        $sq =
                          ";${fs}Query_seq $hsp_seq{$tq} ;${fs}"
                          . "Subject_seq $hsp_seq{$ts}"
                          if $sequence_flg;
                        $prt_funct->();
                        last PRINT;
                    } unless $nogff_flg;

                    #
                    $expanded_flg && do { $prt_funct->(); last PRINT; };

                    #
# ($qm, $sm, $x, $y) = ("$query{$nm}\_\#$n", "$nm\_\#$n", " ", " ");
                    ( $qm, $sm, $x, $y ) = ( "$query{$nm}", "$wnm", " ", " " );
                    $ml = &max( length($qm), length($sm) );
                    ( $a, $b ) = (
                      &fill_right( $qm, $ml, " " ),
                      &fill_right( $sm, $ml, " " )
                    );
                    $prt_funct->();
                };    # PRINT
            };    # for $cnt{$nm}
        };    # do if $cnt{$nm}>0
    };    # foreach
    &prt_foeprg($ht);
    $HT += $ht;
    &header(
        "FROM $. INPUT LINES", "$PT INPUT RECORDS READ",
        "TOTAL HSPs WRITTEN: $HT"
      )
      if $err_flg;
}    # prt_out

sub init_main_vars() {
    @seqlist   = ();
    %prgseq    = ();
    %dbase     = ();
    %query     = ();
    %cnt       = ();
    %desc      = ();
    %sco       = ();
    %hsp_start = ();
    %hsp_end   = ();
    %hsp_seq   = ();
}    # init_main_vars

###################################################
## Main Loop function
###################################################
#
sub check_blastprogname() {
    my ( $str, $pt ) = @_;

    # print STDERR "$_\n";
    ( $prt == 0 && $fragment == 1 ) && do {
        $prt = 1;
        &header( "", "WARNING!!! NO PARAMETERS SECTION FOUND...", "" )
          if $err_flg;
    };
    $prt && do {
        &prt_foeprg($$pt);
        $PT += $$pt;
        $$pt = 0;
        &prt_out;
        &init_main_vars();
        $err_flg && do {
            print STDERR "##\n";
            &header("PARSING ANOTHER BLAST BLOCK");
        };
    };    # $prt
          #
    ( $program, $version ) = split /\s+/o, $str, 2;

    # typeQ/typeS: 0 for proteins - 1 for nucleic acids.
    # Amino Acids vs Amino Acids
    ( $program =~ /^BLASTP$/o ) && ( $scQ = 0, $scS = 0 );

    # Nucleotides vs Nucleotides
    ( $program =~ /^BLASTN$/o ) && ( $scQ = 0, $scS = 0 );

    # Nucleotides vs Amino Acids 
    ( $program =~ /^BLASTX$/o ) && ( $scQ = 1, $scS = 0 );

    # Amino Acids vs Nucleotides
    ( $program =~ /^TBLASTN$/o ) && ( $scQ = 0, $scS = 1 );

    # Nucleotides translated vs Nucleotides translated
    ( $program =~ /^TBLASTX$/o ) && ( $scQ = 1, $scS = 1 );

    #
    $prog_params = "#\n# $program $version\n#";
    $query_name  = $db_name = '';
    $main        = 1;
    $seqflg      = $hsp = $fragment = $param = 0;
}    # check_blastprogname
     #

sub check_seqids() {
    my ( $str, $ifh ) = @_;

    # print STDERR "$_\n";
    defined($program) || do {    # assuming to parse blastn???
        $program = "BLASTN";
        ( $scQ, $scS ) = ( 0, 0 );
        $prog_params = "#\n# $program $version\n#";
        $query_name = $db_name = '';
    };
    ( $seqname, $descr ) = split /\s+/o, $str, 2;
    $seqname =~ s/^\s*>//o;
    $seqname =~ s/\s|:|\|/_/og;
    $seqname .= "_" . ( $index++ );
    $prgseq{$seqname} = "$program ($version)";
    $query{$seqname}  = $query_name;
    ( $db_name =~ /([^\/]+)$/o ) && ( $dbase{$seqname} = $1 );
    push ( @seqlist, $seqname );
    $desc{$seqname} =
      defined($descr) ? join ( ' ', $descr, &get_lines( ' ', $ifh ) ) : '';

    #       $cnt{$seqname} = 0;
    $seqflg = 1;
    $main = $hsp = $fragment = $param = 0;
}    # check_seqids
     #

sub parseblast() {
    my ( $ifh, $file, $pt );
    ( $ifh, $file ) = @_;
    $pt = 0;
    ( $CT, $PT, $HT ) = ( 0, 0, 0 );
    &header( "PARSING BLAST BLOCKS FROM:",
        "  " . ( $file eq '-' ? "STDIN" : $file ) )
      if $err_flg;
    $program = undef;    # default program name
    while (<$ifh>) {
        next if /^\s*$/o;
        &prt_progress( ++$pt );
        my $tmpinput = $_;
        chomp;    # s/\r\n$/\n/;
                  # if your input records finish with "\r\n" (like EMBL).
                  # print STDOUT "$. : $_ \n";
                  # "$." is record number && "$_" is whole record
        CHECK: {

            # Starts with "T?BLAST[PNX]?" ?
            $_ =~ /^\s*T?BLAST[PNX]?/o && do {
                &check_blastprogname( $_, \$pt );
                last CHECK;
            };    # /^\s*T?BLAST[PNX]?/
                  # Starts with ">" ?: sequences.
            $_ =~ /^>/o && do {
                &check_seqids( $_, $ifh );
                last CHECK;
            };    # /^\s*>/
                  # Starts with "Score" ?: HSPs.
            ( $_ =~ /^\s*Score/o && $seqflg ) && do {

                # print STDERR "$_\n";
                ( $score = $_ ) =~ s/^\s*//o;
                $cnt{$seqname}++;
                $sco{ $seqname . $cnt{$seqname} } =
                  join ( '', $score, &get_lines( ', ', $ifh ) );

                # print STDOUT $sco{$seqname.$cnt{$seqname}}."\n";
                $hsp      = 1;
                $fragment = 0;
                last CHECK;
            };    # (/^\s*Score/ && ($seqflg))
                  # Starts with "Query" ?: Fragments.
            ( $_ =~ /^Query:/o && $hsp ) && do {

                # print STDERR "$_\n";
                $fragment = 1;
                last CHECK;
            };    # (/^\s*Query/ && ($hsp))
                  # Parameters Section.
            ( $seqflg && $_ =~ /^\s*(Database|Parameters):/o ) && do {

                # print STDERR "$_\n";
                $prt = $param = 1;
                $seqflg = $hsp = $fragment = 0;
                last CHECK;
            };    # (/^\s*(?:Database|Parameters)/ && ($hsp))
        };   # CHECK Block
             # print STDOUT "$. : MAIN=$main SEQFLG=$seqflg ".
             #              "HSP=$hsp FRAGMENT=$fragment => SEQNAME=$seqname\n";
        LOAD: {
            $fragment && do {    # We are within a fragment.
                my $kk;
                $txt = '';
                $tt  = $cnt{$seqname};
                ( $txt, $kk ) = split /:\s+/o, $_;
                if ( $txt =~ /^Query/o ) {

                    # print STDERR "$_\n";
                    $kk =~ /^(\d+)\s*([^\d]+?)\s*(\d+)\s*$/o && do {
                        ( $ori, $seq, $end ) = ( $1, $2, $3 );
                        ( $hsp_start{ $seqname . "query" . $tt } )
                          || ( $hsp_start{ $seqname . "query" . $tt } = $ori );
                        $hsp_end{ $seqname . "query" . $tt } = $end;
                        $hsp_seq{ $seqname . "query" . $tt } .= $seq;
                    };
                }    # if ($txt =~ /Query/)
                elsif ( $txt =~ /^Su?bje?ct/o ) {

                    # print STDERR "$_\n";
                    $kk =~ /^(\d+)\s*([^\d\s]+)\s*(\d+)\s*$/o && do {
                        ( $ori, $seq, $end ) = ( $1, $2, $3 );
                        ( $hsp_start{ $seqname . "sbjct" . $tt } )
                          || ( $hsp_start{ $seqname . "sbjct" . $tt } = $ori );
                        $hsp_end{ $seqname . "sbjct" . $tt } = $end;
                        $hsp_seq{ $seqname . "sbjct" . $tt } .= $seq;
                        $fragment = 0;
                    };
                }    # elsif ($txt =~ /Sbjct/)
                else { last LOAD; }
            };    # ($fragment)
            $main && do {    # We are within the blast file header.
                /^\s*Query=\s+(.*?)(\s+([^\s]*)?)?$/o && do {

                    # print STDERR "$_\n";
                    ( $query_name = $1 ) =~ s/[\s:\|]/_/g;
                };
                /^\s*Database:\s+(.*?)(\s+([^\s]*)?)?$/o && do {

                    # print STDERR "$_\n";
                    $db_name = $1;
                    while (<$ifh>) {
                        last if /^(?:.*\bletter.*|\s*)$/o;

                        # print STDERR "$_\n";
                        chomp;
                        $_ =~ s/^\s*//o;
                        $_ =~ s/\s*$//o;
                        $db_name .= $_;
                    };    # while getline
                };    # /^\s*Database: +(.*)\s*$/
                last LOAD;
            };    # ($main)
            $param && do {    # We are within the blast file trailer.
                /^\s*Query=\s+(.*?)(\s+([^\s]*)?)?$/o && do {

                    # print STDERR "$_\n";
                    @seqlist = ();
                    ( $query_name = $1 ) =~ s/[\s:\|]/_/g;
                    $main = 1;
                    $seqflg = $hsp = $fragment = $param = 0;
                    last LOAD;
                };
                if (/^\s*[^\[\<\-]/o) {

                    # print STDERR "$_\n";
                    chomp;
                    $_ =~ s/^/\n\# /o;
                    $prog_params = join ( '', $prog_params, $_ );
                }
                else { $param = 0; }
                last LOAD;
            };    # ($param)
        };    # LOAD Block
    };    # while
    $prt && do {
        &prt_foeprg($pt);
        $PT += $pt;
        $pt = 0;
        &prt_out;
        &init_main_vars();
    };    # $prt
          # $prt || &prt_foeprg($pt);
}    # parseblast
