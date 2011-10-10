#!/usr/bin/env perl

=head1 NAME

stamp-split-biounit-models - Splits a biological unit into 1 file per model

=head1 SYNOPSIS

 stamp-split-biounity-models /path/to/biounit


=head1 DESCRIPTION

If you sync the PDB biological asssemblies using the standard file system layout from: 

 rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/biounit/coordinates/

The structures will be gzipped and will contain multiple models per file. As STAMP only processes the first model of each file, this script will split each biological assembly into one file per model. The resulting files are also gzipped and have a file name pattern that identifies the model and the assembly that it came from. E.g. given the structure with accession 1GHD, there are two biological assemblies in the biounit, which are:

 biounit/coordinates/divided/gh/1ghd.pdb1.gz
 biounit/coordinates/divided/gh/1ghd.pdb2.gz

When you run this on C<biounity/coordinates/divided> you will get:

 biounit/coordinates/divided/gh/1ghd-1-1.pdb.gz
 biounit/coordinates/divided/gh/1ghd-2-1.pdb.gz
 biounit/coordinates/divided/gh/1ghd-2-2.pdb.gz

which shows that the first assembly contained only a single model whereas the
second contained two models, which are now in separate files. (The original
files are not modified).

If you want to define a STAMP-readable structure using these models in a
single structure, it is now possible, because they are in separate files,
e.g. in a DOM file my_strucutre.dom:

 /path/to/biounit/coordinates/divided/gh/1ghd-2-1.pdb.gz model1 { ALL }
 /path/to/biounit/coordinates/divided/gh/1ghd-2-2.pdb.gz model2 { ALL }


This requires that you have L<gzip(1)> and L<zcat(1)> installed.

=cut

use strict;
use warnings;
my $dh;
my $fn;
my $fh;
my $fn_out;
my $model;
my $fh_out;
my $dn;
defined($dn = shift @ARGV) or usage();

sub usage {
    my $prog;

    ($prog = __FILE__) =~ s/.*\///;

    die "Usage: $prog pdb_bio_units_directory\n";
}

opendir($dh, $dn) or die("Error: cannot open '$dn' directory.");

while($fn = readdir $dh) {
    # Don't re-read existing model files
    next if $fn =~ /-/;
    # Skip links
    next if -l $fn;
    if($fn =~ /\.gz\Z/) {
	if(!open($fh, "zcat $dn/$fn |")) {
	    warn "Error: cannot open pipe from 'zcat $dn/$fn'";
	    next;
	}
    my $fnctime = (stat("$dn/$fn"))[10];
    my ($pdbid, $assembly) = $fn =~ /^(....)\.pdb(\d+)/;
	while(<$fh>) {
	    if(/^MODEL\s+(\S+)/) {
		$model = $1;
        # Use dash separators, because STAMP strips off extension otherwise
        $fn_out = join('-', $pdbid, $assembly, $model) . '.pdb.gz';
        # Skip unless downloaded file newer
        # (check ctime, since mtime will have been fudged by mirroring)
        my $fn_outctime = (stat("$dn/$fn_out"))[10]; 
        next if ($fn_outctime && $fnctime <= $fn_outctime);
		if(!open($fh_out, "| gzip > $dn/$fn_out")) {
		    warn "Error: cannot open pipe to 'gzip > $dn/$fn_out'";
		    last;
		}
	    }
	    elsif(/^ENDMDL/) {
		close($fh_out) if defined $fh_out;
		$fh_out = undef;
        # Use same file modification time as original
        my @stat = stat("$dn/$fn");
        utime($stat[9], $stat[9], "$dn/$fn_out");
	    }
	    elsif(defined($fh_out)) {
		print $fh_out $_;
	    }
	}
    }
}
closedir($dh);
