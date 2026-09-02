#!/usr/bin/env perl
# ==============================================================================
# Render the flattened HyPR01 mitochondrial GenBank record with OGDraw's
# NATIVE GC-content graph enabled.
#
# OGDraw (GeneMap::Plastome::createMap) supports a `gc_cont => 1` argument that
# draws the standard gray GC-content ring, but the stock `bin/drawgenemap` CLI
# never passes it. This wrapper reproduces exactly what draw_circular_map() in
# drawgenemap does, adding `gc_cont => 1`. It does not modify the OGDraw
# distribution.
#
# Usage (run from the GeneMap-1.1.1 dir with -Ilib, inside ogdraw_env):
#   perl -Ilib render_flattened_with_gc.pl <infile.gb> <outfile.png>
# ==============================================================================
use strict;
use warnings;
use GeneMap::Plastome;

my ($infile, $outfile) = @ARGV;
die "usage: render_flattened_with_gc.pl <infile.gb> <outfile.png>\n"
    unless defined $infile and defined $outfile;

# gc_cont => 1 is a *constructor* argument (it sets the object's _gcGraph flag,
# which createMap then honours by calling _draw_GC_graph).
my $map = GeneMap::Plastome->new(file => $infile, gc_cont => 1);
$map->setTidy(1);
$map->createMap(
    outputfile => $outfile,
    type       => 'png',
    density    => '300x300',
);
print "wrote $outfile (GC graph enabled)\n";
