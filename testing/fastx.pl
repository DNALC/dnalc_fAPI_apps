#!/usr/bin/perl
use strict;
use common::sense;

use lib '/home/smckay/lib';

use AnyEvent;
use IO::File ();
use File::Spec ();
use iPlant::FoundationalAPI::Constants ':all';
use iPlant::FoundationalAPI ();
use Data::Dumper;
use File::Basename;

use constant DEBUG => 0;
use constant iPLANT_USER  => $ENV{USER};
use constant iPLANT_TOKEN => $ENV{TOKEN};

my $cluster = shift || 'stampede';

my $api_instance = iPlant::FoundationalAPI->new(
    debug => DEBUG,
    user  => iPLANT_USER,
    token => shift || iPLANT_TOKEN,
    );

die "Can't auth.." unless $api_instance->auth;
if ($api_instance->token eq kExitError) {
    print STDERR "Can't authenticate!" , $/;								
}

print "Token: ", $api_instance->token, "\n";

my $apps = $api_instance->apps;
my $base = "dnalc-fxtrim-$cluster";
my @apps = $apps->find_by_name($base);

print Dumper [map {$_->id} @apps];

my $cl = pop @apps;

if ($cl) {
    print "Found App ", $cl->name, "\n";
    print STDERR Dumper( $cl ), $/ if DEBUG;
}
else {
    print STDERR  "App [fastx] not found!!", $/;
    exit -1;
}

my $io = $api_instance->io;

my $base_dir = '/' . $api_instance->user;
print "Working in [", $base_dir, "]", $/;

my $job_ep = $api_instance->job;
$job_ep->debug(DEBUG);

my $job_id;

my %params = (
    jobName => 'fx'.int(rand()*1000),
    archive => 1,
    archivePath => "/smckay/API_test/fastx/$cluster\-$$",
    processors => 1,
    requestedTime => '00:30:00',
    softwareName => $cl->name,
    #seq1 => '/smckay/fastq/WT_rep1.fastq',
    seq1 => '/shared/iplant_DNA_subway/sample_data/fastq/small/WT_rep1.fastq',
    #seq1 => '/shared/iplant_DNA_subway/sample_data/fastq/caenorhabditis_elegans/WT_rep1.fastq',
    first => 1,
   'last' => 200,
    quality_threshold => 20,
    min_length => 25,
    min_quality => 20,
    percent_bases => 50
    );

my $job = $job_ep->submit_job($cl, %params), $/;
if ($job != kExitError) {
    $job_id = $job->{data}->{id};
    print STDERR  "JOB_ID: ", $job_id, $/;
}
else {
    print STDERR  "Failed to submit job..", $/;
}

unless ($job_id) {
    die "Job not submitted..\n", Dumper $job;
}

print STDERR  "Polling for job status..", $/;

my $i = 20;

while ($i) {
    my $st = $job_ep->job_details($job_id);
    my $stat = $job_id. "\t". $st->{status}. "\t" . `date`;
    print $stat;
    last if $stat =~ /FINISHED/;
    $i--;
    sleep 30;
}

