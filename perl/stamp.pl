#!/usr/bin/perl
use strict;
use warnings;
use Carp;
use English;
use FindBin;
use Getopt::Long;
use lib "$FindBin::Bin/../lib/perl";
use PDB::Entry;
use STAMP::Domain;
use STAMP::Tools;
use constant PDB_ID_LENGTH => 4;

if ( !$ENV{STAMPDIR} ) {
	croak "error: environment variable STAMPDIR is not defined\n";
}

my $pdbcodes;
my $domain_file;
my $list;
my $help;
my $align;
my $prefix = 'stamp';

GetOptions(
    'align'         => \$align,
	'domain-file=s' => \$domain_file,
    'help'          => \$help,
	'list'          => \$list,
	'pdb=s'         => \$pdbcodes,
);

my @descriptors;
if ($pdbcodes) {
	my @pdbcodes = split /,/, $pdbcodes;

	for my $pdbid (@pdbcodes) {
		my $id = substr $pdbid, 0, PDB_ID_LENGTH;
		my $pdbfile = getfile($pdbid)
			or croak "cannot find PDB file for $pdbid";

		my $pdb = PDB::Entry->new($pdbfile);
		for my $chain ( sort keys %{ $pdb->chains } ) {
			my $desc      = "$id $id$chain { CHAIN $chain }\n";
			my $domain_ok = 1;
			for my $residue ( @{ $pdb->chains->{$chain} } ) {
				if ( $residue->is_amino && !$residue->atom('CA') ) {
					my $to_str = $residue->to_string;
					warn <<"WARNING";
% WARNING Residue $to_str in chain $chain of $pdbfile is missing its mainchain CA atom. 
% STAMP cannot use this chain. 
WARNING
					$domain_ok = 0;
				}
			}
			if ($domain_ok) {
				push @descriptors, $desc;
			}
		}
	}
}

if ($domain_file) {
	open my $IN, '<', $domain_file or croak "$!: $domain_file";
	while ( my $line = <$IN> ) {
		chomp $line;
		next if $line =~ /^%/;
		my $domain  = STAMP::Domain->new($line);
		my $pdbfile = getfile( $domain->pdb );
		if ( !$pdbfile ) {
			croak "Cannot find PDB file for domain descriptor $line";
		}
		my $desc;
		if (eval { $desc = STAMP::Tools::descriptor_ok( $domain, $pdbfile ); }){
			push @descriptors, "$desc\n";
		}
		else {
			my $error = ( split /\n/, $EVAL_ERROR )[0];
			die <<"ERROR";
% Error in descriptor: $line
$error
ERROR
		}
	}
	close $IN;
}

if ($list) {
	print @descriptors;
}

if ($align){
    do_align(\@descriptors);
}

sub do_align {
    my $descriptors = shift;
    my ($first,@others) = @$descriptors;
    my $prefix = $first;
    $prefix =~ tr/[ {}]/_/s;
    chomp $prefix;
    mkdir $prefix;
    chdir $prefix;
    my $query_file = "$prefix.query";
    open my $QFH,'>',$query_file or die "$!: $query_file";
    print {$QFH} $first;
    close $QFH;
    my $db_file = "$prefix.db";
    open my $DFH, '>', $db_file or die "$!: $db_file";
    print {$DFH} @others;
    close $DFH;

    system "stamp -s -d $db_file -f $query_file -prefix $prefix";
    system "sorttrans -f $prefix.scan > $prefix.out";
    system "stamp -f $prefix.scan -prefix $prefix";
    # Suffix for final STAMP MSA file
    my $suffix = @others - 2;
    my $bloc_file = "$prefix.$suffix";
    system "aconvert -in b -out f <  $bloc_file > $prefix.fa";
}

sub getfile {
	my ( $pdbid, $type ) = @_;
	if (! defined $type ) {
		$type = 'pdb';
	}
	my $dirfile = "$ENV{STAMPDIR}/$type.directories";

	my $file;
	open my $FH, '<', $dirfile or croak "$!: $dirfile";
	while (<$FH>) {
		my ( $path, $prefix, $suffix ) = map { $_ eq q{_} ? q{} : $_ } split;
		if ($path) {
			$path .= q{/};
		}
		$file = "$path$prefix$pdbid$suffix";
		if ( -e $file ) {
			last;
		}
		undef $file;
	}
	close $FH;
	return $file;
}
