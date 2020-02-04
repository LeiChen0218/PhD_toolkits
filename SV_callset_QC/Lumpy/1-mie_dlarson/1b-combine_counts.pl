#!/usr/bin/perl

use strict;
use warnings;
my %totals;
while(<>) {
    chomp;
    my ($trio, $type, $filter, $msq, $mie_class, $count) = split "\t";
    my $msq_bin = int($msq/10)*10;
    my $key = join("\t", $type, $filter, $msq_bin);
    $totals{$key}{$mie_class} += $count;
}

print join("\t", qw(SVType Filter MSQBin Number MIE NoMIE MIERate)), "\n";
for my $key (keys %totals) {
    my $total = 0;
    for my $class qw( MIE NONE) {
        if (exists($totals{$key}{$class})) {
            $total += $totals{$key}{$class};
        }
        else {
            $totals{$key}{$class} = 0;
        }
    }
    my $mie_rate = 'NA';
    if ($total > 0) {
        $mie_rate = $totals{$key}{'MIE'} / ($totals{$key}{'NONE'} + $totals{$key}{'MIE'});
    }

    print join("\t", $key, $total, $totals{$key}{'MIE'}, $totals{$key}{'NONE'}, $mie_rate), "\n";
}
