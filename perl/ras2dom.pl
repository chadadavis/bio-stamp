#!/usr/local/bin/perl -w

$code = "unk";

if(defined($ARGV[0])) { $code = $ARGV[0] }

$n=0;
$ndom=1;
while(<stdin>) {
   if(/^Atom:/) {
       chop;
       @t = split(/ +/);
       $chain = "_";
       $ins = "_";
       if(/Chain/) { $chain = $t[7] }
       if($n==0) { 
          $start = $chain . " " . $t[5] . " " . $ins;
          undef($end);
#       } elsif((($n+1)%2)==0) { 
       } else { 
          $end = $chain . " " . $t[5] . " " . $ins;
          printf("UNK %s_%-2d { %s to %s }\n",$code,$ndom,$start,$end);
          $ndom++;
          $start = $end;
       }
#       print "n= ",$n," ",$_,"\n";
       $n++;
       
    }
}

