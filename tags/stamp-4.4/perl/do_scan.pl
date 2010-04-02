#!/usr/local/bin/perl -w

$|=1;

$domfile = $ARGV[0];

$id = $domfile;
$id =~ s/\.domain//;
$id =~ s/\.dom//;

if((!defined($domfile)) || (!-e $domfile)) {
   die "You must specify a domain file name\n";
}

$options = " -n 2 -opd -secscreen F ";
$clean_value = 3;
$probefile = "/tmp/probe.".$$.".dom";
$dbasefile = "/tmp/dbase.".$$.".dom";

open(DOM,$domfile)     || die "Error reading file $domfile\n";
open(PROBE,">$probefile") || die "Error writing to file $probefile\n";
open(DBASE,">$dbasefile") || die "Error reading file $dbasefile\n";

$ndom = 0;
while(<DOM>) {
   chop;
   if((!/^[%\#]/) && (/{/)) {
      if($ndom ==0) { 
	 print PROBE $_; 
	 if(!/}/) { print PROBE "}"; } 
	 print "PROBE will be ",(split(/ +/))[1],"\n";
	 print PROBE "\n";
	 close(PROBE);
      } else {
	 print DBASE $_; 
	 if(!/}/) { print DBASE "}"; } 
	 print " compared to ",(split(/ +/))[1],"\n";
	 print DBASE "\n";
      }
      $ndom++;
   }
}
close(DBASE);
close(DOM);

$transfile = $id . ".trans";
open(TRANS,">$transfile") || die "Error writing to file $transfile\n";

$com = "stamp -l ".$probefile." -s -d ".$dbasefile . " -prefix /tmp/stamp.".$$ . $options;
open(IN,"$com|") || die "Error running $com\n";
while(<IN>) {
   if(/^Scan/) { print TRANS "%",$_; }
   print;
}
close(IN);


$scanfile = "/tmp/stamp.".$$.".scan";
open(IN,$scanfile) || die "Error reading $scanfile\n";
while(<IN>) {
   if(!/^[%#]/) { 
      s/_[0-9][0-9 ]*/ /;
      print TRANS; 
   }
} 
close(TRANS);

$com = "stamp -l ".$transfile. " -n 2 -prefix /tmp/stamp.";
open(IN,"$com|") || die "Error running $com\n";
while(<IN>) {
   print;
   if(/^ See file .* for the alignment and transformations/) {
       s/^ //;
       $lastfile = (split(/ +/))[2];
   }
}
close(IN);

$alignfile = $id . ".stamp";
$com = "stamp_clean ".$lastfile." ".$clean_value." > ".$alignfile;
system($com);

print "Transformations derived from scanning should be in $transfile\n";
print "Final alignment/superimposition in $alignfile\n";

unlink $probefile;
unlink $dbasefile;

$com = "/bin/rm /tmp/stamp.".$$."*";
system($com);


