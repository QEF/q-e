#!/usr/bin/perl -w

use strict;

{
    my $maxdepth = 2; # default
    if ($#ARGV > 0 && $ARGV[0] eq "-d")
    {
	$maxdepth = $ARGV[1];
	if ($maxdepth !~ /^\d+$/)
	{
	    print STDERR "error: non-numeric maxdepth\n";
	    print STDERR "usage: $0 [-d maxdepth] [targets]\n";
	    exit 1;
	}
	shift @ARGV;
	shift @ARGV;
    }

    # $basedir is directory where this script is
    my $basedir = $0;
    $basedir =~ s/(.*)\/.*/$1/;
    my @sources = split(/ /, `echo $basedir/*/*.f90`);

    # grab program, function and subroutine declarations
    my (%place, %fname, %pname);
    foreach my $file (@sources)
    {
	open(IN, "$file");
	while (<IN>)
	{
	    $_ = "\L$_"; # cast everything to lowercase
	    if (/^[^!'""']*\bfunction\s+(\w+)/o && ! /^\s*end\s+function\b/o)
	    {
		$fname{$1} = 1;
		push_place(\%place, $1, $file);
	    }
	    elsif (/^\s*program\s+(\w+)/o)
	    {
		$pname{$1} = 1;
		push_place(\%place, $1, $file);
	    }
	    elsif (/^\s*(?:(?:pure|recursive)\s+)?subroutine\s+(\w+)/o)
	    {
		push_place(\%place, $1, $file);
	    }
	}
	close(IN);
    }
    my @names = sort keys %place;
    my @functions = sort keys %fname;

    # if no arguments are specified, stat all programs
    my @targets = @ARGV;
    if ($#targets < 0)
    {
	@targets = sort keys %pname;
    }

    my %cache;
    foreach my $name (@targets)
    {
	stat_name($name, \%place, \@functions, \%cache, 0, "", $maxdepth);
    }
}

sub push_place
{
    my ($place, $name, $file) = @_;
    if (defined $$place{$name})
    {
	$$place{$name} = "$$place{$name} $file";
    }
    else
    {
	$$place{$name} = "$file";
    }
}

sub stat_name
{
    my ($name, $place, $functions, $cache, $depth, $indent, $maxdepth) = @_;
    print "$indent$name\n";

    if ($depth >= $maxdepth || ! defined $$place{$name})
    {
	return;
    }

    if (! defined $$cache{$name})
    {
	my %cname;
	my @files = split(/ /, $$place{$name});
	foreach my $file (@files)
	{
	    my $current = "";
	    open(IN, $file);
	    while (<IN>)
	    {
		$_ = "\L$_";
		if (/^\s*program\s+(\w+)/o)
		{
		    $current = "$1";
		}
		elsif (/^\s*(?:(?:pure|recursive)\s+)?subroutine\s+(\w+)/o)
		{
		    $current = "$1";
		}
		elsif (/^[^!'""']*\bfunction\s+(\w+)/o)
		{
		    $current = "$1";
		}
		# here we are inside the relevant program/subroutine/function
		elsif ($current eq $name)
		{
		    # subroutine calls
		    if (/^\s*call\s+(\w+)/o)
		    {
			$cname{$1} = 1;
		    }
		    # function calls
		    foreach my $fun (@$functions)
		    {
			if (/^[^!'""']*\b$fun\b/)
			{
			    $cname{$fun} = 1;
			}
		    }
		}
	    }
	    close(IN);
	}
	my @calls = sort keys %cname;
	$$cache{$name} = \@calls;
    }

    foreach my $call (@{$$cache{$name}})
    {
	if ($call ne $name)
	{
	    stat_name($call, $place, $functions, $cache,
		      $depth+1, "  $indent", $maxdepth);
	}
    }
}
