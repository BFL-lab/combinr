#!/usr/bin/env perl
# Cross-check harness: run PASA's real CDNA::Alternative_splice_comparer on the
# isoform structures emitted by combinr (parsed from its isoform GFF3), and print
# canonical events to compare against combinr's own event report.
#
# Usage: pasa_altsplice.pl <PASA_PerlLib_dir> <combinr_isoforms.gff3>

use strict;
use warnings;

my ($perllib, $gff3) = @ARGV;
die "usage: $0 <PASA/PerlLib> <isoforms.gff3>\n" unless $perllib && $gff3;

# Add the PASA PerlLib to @INC at runtime, then load its modules.
require lib;
lib->import($perllib);
require Gene_obj;
require CDNA::Alternative_splice_comparer;

# ---- parse combinr GFF3: gene -> ordered list of isoforms {id, strand, exons} ----
my %gene_order;          # gene_id -> [transcript_id,...] in document order
my %tx;                  # transcript_id -> { gene, strand, exons=>[[l,r],...] }
open(my $fh, '<', $gff3) or die "cannot open $gff3: $!";
while (<$fh>) {
    next if /^#/;
    chomp;
    my @c = split /\t/;
    next unless @c >= 9;
    my ($type, $lend, $rend, $strand, $attr) = ($c[2], $c[3], $c[4], $c[6], $c[8]);
    if ($type eq 'mRNA') {
        my ($id)     = $attr =~ /ID=([^;]+)/;
        my ($parent) = $attr =~ /Parent=([^;]+)/;
        push @{$gene_order{$parent}}, $id;
        $tx{$id} = { gene => $parent, strand => $strand, exons => [] };
    } elsif ($type eq 'exon') {
        my ($parent) = $attr =~ /Parent=([^;]+)/;
        push @{$tx{$parent}{exons}}, [$lend, $rend];
    }
}
close $fh;

sub build_gene {
    my ($id) = @_;
    my $t = $tx{$id};
    my $g = new Gene_obj();
    $g->build_gene_obj_exons_n_cds_range($t->{exons}, 0, 0, $t->{strand});
    return $g;
}

my $cmp = new CDNA::Alternative_splice_comparer();
my @out;
sub emit {
    my ($kind, $coords, $a, $b, $detail) = @_;
    push @out, join("\t", $kind, $coords, $a, $b, $detail // '');
}

sub coordstr {
    # sort coordpairs by lend, join as "l-r,l-r"
    my @segs = sort { $a->[0] <=> $b->[0] } @_;
    return join(",", map { "$_->[0]-$_->[1]" } @segs);
}

for my $gene (sort keys %gene_order) {
    my @isos = @{$gene_order{$gene}};
    for (my $i = 0; $i < @isos; $i++) {
        for (my $j = $i + 1; $j < @isos; $j++) {
            my ($id1, $id2) = ($isos[$i], $isos[$j]);
            my $g1 = build_gene($id1);
            my $g2 = build_gene($id2);

            # retained introns: both directions
            for my $p ([$g1, $g2, $id1, $id2], [$g2, $g1, $id2, $id1]) {
                my ($ga, $gb, $ia, $ib) = @$p;
                for my $in ($cmp->find_unspliced_introns($ga, $gb)) {
                    emit("retained_intron", "$in->[0]-$in->[1]", $ia, $ib, "");
                }
            }

            # conventional alt donor/acceptor: once (a,b)
            my %conv = $cmp->find_conventional_alt_splice_isoforms($g1, $g2);
            for my $acc (@{$conv{acceptors} || []}) {
                my ($x, $y) = sort { $a <=> $b } ($acc->{gene1}, $acc->{gene2});
                emit("alt_acceptor", "$x-$y", $id1, $id2, "");
            }
            for my $don (@{$conv{donors} || []}) {
                my ($x, $y) = sort { $a <=> $b } ($don->{gene1}, $don->{gene2});
                emit("alt_donor", "$x-$y", $id1, $id2, "");
            }

            # starts/ends within introns: both directions
            for my $p ([$g1, $g2, $id1, $id2], [$g2, $g1, $id2, $id1]) {
                my ($ga, $gb, $ia, $ib) = @$p;
                my %se = $cmp->find_starts_and_ends_within_introns($ga, $gb);
                emit("start_within_intron", "$se{start}-$se{start}", $ia, $ib, "") if exists $se{start};
                emit("end_within_intron",   "$se{end}-$se{end}",     $ia, $ib, "") if exists $se{end};
            }

            # exon skipping: both directions
            for my $p ([$g1, $g2, $id1, $id2], [$g2, $g1, $id2, $id1]) {
                my ($ga, $gb, $ia, $ib) = @$p;
                for my $event ($cmp->find_exon_skipping_events($ga, $gb)) {
                    emit("exon_skip", coordstr(@$event), $ia, $ib, "");
                }
            }

            # alternate terminal exons: both directions
            for my $p ([$g1, $g2, $id1, $id2], [$g2, $g1, $id2, $id1]) {
                my ($ga, $gb, $ia, $ib) = @$p;
                for my $alt ($cmp->find_alternate_exons($ga, $gb)) {
                    my $side = $alt->{type} eq 'lend' ? 'front' : 'back';
                    my ($rl, $rr) = @{$alt->{coords}};
                    emit("alternate_exon", "$rl-$rr", $ia, $ib,
                         "side=$side;num_exons=$alt->{num_exons}");
                }
            }
        }
    }
}

# de-duplicate identical rows (combinr dedups too) and print sorted
my %seen;
print "$_\n" for grep { !$seen{$_}++ } sort @out;
