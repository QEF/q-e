#!/usr/bin/perl -w

use strict;

{
    # $basedir is directory where this script is
    my $basedir = $0;
    $basedir =~ s/(.*)\/.*/$1/;
    my @sources = split(/ /, `echo $basedir/*/*.f90`);

    # grab program, function and subroutine declarations
    my (%place, %fname, %pname, %sname);
    foreach my $file (@sources)
    {
	open(IN, "$file");
	while (<IN>)
	{
	    $_ = "\L$_"; # cast everything to lowercase
	    if (/^[^!'""']*\bfunction\s+(\w+)/o && ! /^\s*end\s+function\b/o)
	    {
		$fname{$1} = 1;
		insert_place(\%place, $1, $file);
	    }
	    elsif (/^\s*program\s+(\w+)/o)
	    {
		$pname{$1} = 1;
		insert_place(\%place, $1, $file);
	    }
	    elsif (/^\s*(?:(?:pure|recursive)\s+)?subroutine\s+(\w+)/o)
	    {
		$sname{$1} = 1;
		insert_place(\%place, $1, $file);
	    }
	}
	close(IN);
    }
    my @targets = sort keys %place;
    my @programs = sort keys %pname;
    my @functions = sort keys %fname;

    # html preamble
    print "<html>\n";
    print "<body>\n";
    print "\n";

    # list of programs
    print "<dl>\n";
    print "  <dt>list of programs:</dt>\n";
    print "  <dd><p>\n";
    foreach my $program (@programs)
    {
	print "    <a href=\"#$program\">$program</a>\n";
    }
    print "  </p></dd>\n";
    print "</dl>\n";
    print "\n";

    # list of all routines
    print "<dl>\n";
    foreach my $name (@targets)
    {
	print "  <dt><a name=\"$name\">";
	if (defined $pname{$name})
	{
	    print "program ";
	}
	elsif (defined $sname{$name})
	{
	    print "subroutine ";
	}
	elsif (defined $fname{$name})
	{
	    print "function ";
	}
	print "$name</a></dt>\n";

	my %cname;
	my @files = split(/ /, $place{$name});
	foreach my $file (@files)
	{
	    print "  <dd><p>defined in file: $file<br>\n";
	    print "    calls:\n";

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
		    foreach my $fun (@functions)
		    {
			if ($fun ne $name && /^[^!'""']*\b$fun\b/)
			{
			    $cname{$fun} = 1;
			}
		    }
		}
	    }
	    close(IN);
	    my @calls = sort keys %cname;
	    foreach my $call (@calls)
	    {
		print "    <a href=\"#$call\">$call</a>\n";
	    }
	    print "  </p></dd>\n";
	}
    }
    print "</dl>\n";

    # html postamble
    print "\n";
    print "</body>\n";
    print "</html>\n";
}

sub insert_place
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
