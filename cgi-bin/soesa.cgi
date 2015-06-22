#!/usr/bin/perl -w
# This is a toy cgi script that can be used to provide a web interface
#   to structure evaluation with soesa.

# Include CGI module:

use CGI qw/:standard/;

$CGI::POST_MAX=1024 * 100;  # max 100K posts

# The directory where work files will be written:

$rootdir = "/tmp/pdf/";
$evalfile = $rootdir."eval.dat";
$exptdb = $rootdir . "expts";            # File for holding the expt entries

$soesahome = "[Place soesa root directory here]";
$soesaopts = "-aacfg ".$soesahome."/dat/aa.cfg -hashfile ".$soesahome."/dat/pdf.hash -datafile ".$soesahome."/dat/pdf.data -evalpdb -out ".$evalfile;
# This subroutine strips off the leading path info from a filename:

sub basename {
    my $agent;
    my @fparse;
    my @fparserev;
    my $fn;

    my $fname = $_[0];

    my $msg = "";

# Action depends on the operating system.  UNIX uses '/', Windows '\':
#   This uses the cryptic regular expression syntax, which resembles
#   vi search and replace syntax, sed syntax, etc.  Seems like you can do
#   most anything with it but easily understand it once it's written:

    $agent=user_agent();
    if ($agent =~ /Win/) {
	$fname =~ s/\\/\//g;
    }
    if (($agent =~ /IRIX/) || ($agent =~ /Linux/) || ($agent =~ /Win/)) {
	@fparse = split /\//,$fname;
	@fparserev = reverse @fparse;
	$fn = $fparserev[0];
    } else {
	print "Cannot recognize browser type.\n",p;
	die;
    }
    return $fn;
}



# Begin the html document:

print header;
if (!param() || param('Refresh')) {
    print start_html('SOESA');
}
if (param('Evaluate')) {
    print start_html('SOESA: Evaluate a structure');
}
if (param('Add')) {
    print start_html('SOESA: Add a structure entry');
}
if (param('Delete')) {
    print start_html('SOESA: Delete structure entries');
}

# Code to evaluate the protein structure:


if (param('Evaluate')) { 
    my $entrynum = param('entryselect');
    
# The filename that holds the entries is $exptdb (defined above):

    if (stat($exptdb)) {

# In perl, read/write/append has syntax like redirection in unix.
#   This opens the file for reading:

	open (FOO, "<$exptdb");

# This is a trick for breaking a text file into lines:

	@lines = <FOO>;
	close FOO;

# If the 'All' option was selected, give an error msg:

	if ($entrynum == 0) {
	    print "'All' is not a valid option for structure evaluation.",br;
	}

# Otherwise evaluate the selected pdb:

	else {
	    $soesacmd = $soesahome."/bin/soesa ".$soesaopts." -pdb ".$lines[$entrynum-1]."&";
	    print $soesacmd,br;
	    `$soesacmd`;
	}

    }
}

# The following sections (Delete, Add) relate to maintenance of a list of
#   experiment entries that figure in calculation of the bsf.

# Delete an entry, or all entries from the list.  Also delete associated
#   uploaded data files:

if (param('Delete')) {
    my $entrynum = param('entryselect');

# The filename that holds the entries is $exptdb (defined above):

    if (stat($exptdb)) {

# In perl, read/write/append has syntax like redirection in unix.
#   This opens the file for reading:

	open (FOO, "<$exptdb");

# This is a trick for breaking a text file into lines:

	@lines = <FOO>;
	close FOO;

# If all entries are deleted, also delete the entry file:

	if (($entrynum == 0) || (@lines == 1)) {
	    foreach $line (@lines) {
		unlink $line;
	    }
	    unlink $exptdb;
	    print "All entries deleted.",br;
	}
# Otherwise just remove the entry from the list and delete the data file:

	else {
	    $linenum = $entrynum - 1;
	    $line = $lines[$linenum];
	    unlink $line;
	    print "Entry ",basename($line)," deleted.";

# The splice() function is a handy perl subroutine for removing entries
#   from a list.  This command removes a single entry at position $linenum:

	    splice(@lines,$linenum,1);
	    unlink $exptdb;
	    open (FOO,">$exptdb");
	    foreach $line (@lines) {
		print FOO $line;
	    }
	}
    }
}

# Subroutine to add an entry to the list:

if (param('Add')) {
    my $fname;

# A file must be uploaded.  In the CGI syntax, once the filename is obtained,
#   it can be opened and read from as if it were local to the machine:

    if (!($fname = param('filename'))) {
	print "Please supply a file name.\n",p;
    } else {

	$fn = basename($fname);

# Build the real filename:

	$fn = $rootdir . $fn;

# Upload the file:

	if (!(stat($fn))) {
	    open (FOO, ">$fn") or 
		die "Can't open file.";
	    while (read($fname,$data,1024)) {print FOO $data;}
	    close(FOO);
	}
	else {
	    print "File already exists.",br;
	}

	$msg = "File successfully added.\n";

    }

# Append the info for this entry to the entry file:

    open (FOO, ">>$exptdb");
    
    print FOO $fn,"\n";
    
    close (FOO);

}

# Display and set up menu options:

print start_form(-enctype=>'multipart/form-data');
if (stat($evalfile)) {
    print "<a href=http://milford.mst.lanl.gov:8080/pdfresults/",basename($evalfile),"> Download</a> the evaluation results.",br;
}
print " PDB file:",filefield('filename'),br; 
print submit('Add')," the above entry to the list.",br;
print hr;

if (stat($exptdb)) {
    my $displine;
    my $linenum = 1;
    my @radiolist = ();
    my %radiolabel = ();
    print submit('Evaluate')," the selection below:",br;
    print submit('Delete')," the selection below:",br;
    open (FOO, "<$exptdb");
    @lines = <FOO>;
    push(@radiolist,0);
    $radiolabel{0} = "All entries";
    foreach $line (@lines) {
	push(@radiolist,$linenum);
	$radiolabel{$linenum} = basename($line);
	$linenum++;
    }
    print radio_group('entryselect',\@radiolist,'','true',\%radiolabel);
    close FOO;
}

print hr;

print submit('Refresh')," this page.",br;

print end_form;

print hr;

print end_html;
