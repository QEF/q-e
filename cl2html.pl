#!/usr/bin/perl
# cl2html.pl - Print CVS activity in HTML format based on cvs2cl.pl XML output

# Copyright (C) 2004 Simon Josefsson
# Copyright (C) 2003 Anderson Lizardo
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA


# This program is based on cl2html.pl, retrieved from:
#
# http://www.linuxfromscratch.org/~jeroen/cl2html.pl
#
# but rewritten to produce output that is similar to the output on:
#
# http://curl.haxx.se/auto/cvshistory.html
#
# For updates to this project, see:
#
# http://josefsson.org/cl2html/

use strict;
use warnings;

use XML::Parser;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(strftime);

my $help = 0;
my $version = 0;
my $entries_limit = 50;
my $viewcvs_url = "";
my $viewcvs_extra = "";

GetOptions(
    "help"		=> \$help,
    "version"		=> \$version,
    "viewcvs-url:s"	=> \$viewcvs_url,
    "viewcvs-extra:s"	=> \$viewcvs_extra,
    "entries:i"		=> \$entries_limit
) or pod2usage(1);

pod2usage(1) if $help;

# Handle --verison.
if ($version) {
    my ($rev) = '$Revision: 1.1 $';
    $rev = $1 if $rev =~ m/Revision: (\S+) /;
    print "cl2html $rev\n";
    print "\n";
    print "Copyright (C) 2004 Simon Josefsson\n";
    print "Copyright (C) 2003 Anderson Lizardo\n";
    print "This is free software; see the source for copying conditions.  There is NO\n";
    print "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n";
    exit 0;
}

$viewcvs_extra = "&$viewcvs_extra" if $viewcvs_extra;

my $buffer = "";     # Current text in buffer.
my $entry_count = 0; # Number of entries printed.
my $last_date = "";  # Last date YYYY-MM-DD.
my $last_file = "";  # Last filename, used for viewcvs revision URLs.

# Convert ISO 8601 date (yyyy-mm-ddThh:mm:ssZ) to the specified format
sub isodate2any {
    my ($date, $format) = @_;
    if ($date =~ /(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z/) {
	return strftime($format, $6, $5, $4, $3, $2 - 1, $1 - 1900);
    } else {
	return undef;
    }
}

sub handle_StartTag {
    my (undef, $tag) = @_;

    if ($tag eq "changelog") {
	print "<table border=0>\n";
    }

    $buffer = "";
}

sub handle_EndTag {
    my (undef, $tag) = @_;

    if ($tag eq "tagisodate") {
	print "<tr bgcolor=lightgreen>\n";
	print "<td>" . isodate2any($buffer, '%H:%M') . "\n";
    }
    elsif ($tag eq "tagdatetag") {
	print "<td><td>\n";
	print "<td>tag <b>" . $buffer . "</b> added\n";
	print "</tr>\n";
    }
    elsif ($tag eq "isoDate") {
	if (isodate2any($buffer, '%Y-%m-%d') ne $last_date) {
	    $last_date = isodate2any ($buffer, '%Y-%m-%d');
	    print "<tr>\n";
	    print "<td colspan=4 bgcolor=lightblue><b>$last_date</b>\n";
	}
	if ($entry_count % 2) {
	    print "<tr>\n";
	} else {
	    print "<tr bgcolor=lightgrey>\n";
	}
	print "<td>" . isodate2any($buffer, '%H:%M') . "\n";
    }
    elsif ($tag eq "author") {
	print "<td>$buffer\n";
	print "<td>\n";
    }
    elsif ($tag eq "name") {
	$last_file = $buffer;
    }
    elsif ($tag eq "revision") {
	if ($viewcvs_url) {
	    print "<a href=\"$viewcvs_url/$last_file?"
		."rev=$buffer&view=auto$viewcvs_extra\">"
		. "$last_file</a>\n";
	} else {
	    print "$last_file\n";
	}
	my ($prev, $lastrev);
	if ($buffer =~ /([0-9.]+\.)([0-9]+)/) {
	    $prev = $2 - 1;
	    $lastrev = $1.$prev if $prev != 0;
	}
	if ($viewcvs_url && $lastrev) {
	    print "<a href=\"$viewcvs_url/$last_file?"
		."r1=$lastrev&r2=$buffer$viewcvs_extra\">"
		. "$buffer</a><br>\n";
	} else {
	    print "$buffer<br>\n";
	}
    }
    elsif ($tag eq "msg") {
	$buffer =~ s/\n\n/<br>/gs;
	print "<td>$buffer\n";
    }
    elsif ($tag eq "entry") {
	print "</tr>\n\n";
	if ($entries_limit && ++$entry_count == $entries_limit) {
	    print "</table>\n";
	    print "<p>This page was created on " . (scalar localtime) .
		" using <a href=\"http://josefsson.org/cl2html/\">" .
		"cl2html</a> written by Simon Josefsson.\n";
	    exit 0;
	}
    }
    elsif ($tag eq "changelog") {
	print "</table>\n";
	print "<p>This page was created on " . (scalar localtime) .
	    " using <a href=\"http://josefsson.org/cl2html/\">" .
	    "cl2html</a> written by Simon Josefsson.\n";
    }
}

sub handle_Text {
    my ($expat, $text) = @_;

    # Encode "special" entities
    $text =~ s/\&/\&amp;/g;
    $text =~ s/</\&lt;/g;
    $text =~ s/>/\&gt;/g;
    #$text =~ s/\"/\&quot;/g;
    #$text =~ s/\'/\&apos;/g;

    # Add current text to the buffer
    $buffer .= $text;
}

my $parser = new XML::Parser(
    Handlers => {
	Start => \&handle_StartTag,
	End => \&handle_EndTag,
	Char => \&handle_Text,
    },
);

$parser->parse(\*STDIN);

__END__

=head1 NAME

cl2html - Print CVS activity in HTML format based on cvs2cl.pl XML output

=head1 SYNOPSIS

cl2html  [--help] [--entries NUMBER] [--viewcvs-url BASEURL]

--entries NUMBER     Print at most NUMBER of entries, use 0 for all.
--viewcvs-url URL    Create links to viewcvs.cgi output for each file.
--viewcvs-extra URL  Additional CGI parameters to add to URLs.

--help               Show brief help message.
--version            Display version number.

=head1 DESCRIPTION

B<cl2html> converts the XML outputted by cvs2cl.pl's C<--xml>
option to HTML code.  The data is read on standard input.

=head1 OPTIONS

=over

=item B<--entries NUMBER>

Specify the number of log entries to print.  Use 0 to print all log
entries.  By default, 50 entries are printed.

=item B<--viewcvs-url URL>

Specify location of viewcvs.cgi script that serve CVS content.  Used
to create links to the files, and revisions, in the output.  For
example, you could use:

    --viewcvs-url http://josefsson.org/cgi-bin/viewcvs.cgi/libidn/

=item B<--viewcvs-extra PARAM>

Specify optional extra viewcvs.cgi parameters to be added to
viewcvs.cgi URLs.  This can be used when your real repository looks,
for example, like:

  http://josefsson.org/cgi-bin/viewcvs.cgi/?root=gnupg-mirror

To get this to work, you would use:

    --viewcvs-root http://josefsson.org/cgi-bin/viewcvs.cgi/
    --viewcvs-extra root=gnupg-mirror

=item B<--help>

Print a brief help message and exits.

=item B<--version>

Print the version number and exits.

=back

=head1 AUTHOR

Copyright (C) 2004 Simon Josefsson
Copyright (C) 2003 Anderson Lizardo

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

=cut
