//! Evidence-derived splice / start / stop site sets.
//!
//! EVM scans the genome for canonical GT-AG donors/acceptors and ATG/stop codons to
//! gate candidate exons (`populate_splice_sites`, `populate_starts_and_stops`). The
//! consensus path deliberately does NOT (see `avoid-canonical-splice-bias`): the only
//! admissible splice boundaries are those *observed in the input evidence* (the union
//! of intron boundaries across all chains), and candidate starts/stops are the input
//! predictions' CDS bounds (plus the configurable genetic-code stops found inside
//! candidate exons at frame-enumeration time — those are handled in
//! [`super::exon::determine_good_phases`], not stored here).

use std::collections::HashSet;

/// Sets of evidence-supported sites, in forward genomic coordinates.
#[derive(Default, Debug)]
pub struct SiteSets {
    /// Intron 5' boundary (first base of an intron = exon end3 + 1).
    pub donors: HashSet<i64>,
    /// Intron 3' boundary (last base of an intron = next exon lend - 1).
    pub acceptors: HashSet<i64>,
    /// Prediction CDS start positions (5' base of a start codon).
    pub starts: HashSet<i64>,
    /// Prediction CDS stop positions (3' base of the CDS).
    pub stops: HashSet<i64>,
}

impl SiteSets {
    /// Record the boundaries of an evidence intron `[intron_lend, intron_rend]`.
    pub fn add_intron(&mut self, intron_lend: i64, intron_rend: i64) {
        self.donors.insert(intron_lend);
        self.acceptors.insert(intron_rend);
    }

    pub fn add_start(&mut self, g: i64) {
        self.starts.insert(g);
    }

    pub fn add_stop(&mut self, g: i64) {
        self.stops.insert(g);
    }

    pub fn is_donor(&self, g: i64) -> bool {
        self.donors.contains(&g)
    }

    pub fn is_acceptor(&self, g: i64) -> bool {
        self.acceptors.contains(&g)
    }

    pub fn is_start(&self, g: i64) -> bool {
        self.starts.contains(&g)
    }

    pub fn is_stop(&self, g: i64) -> bool {
        self.stops.contains(&g)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intron_boundaries_become_donor_acceptor() {
        let mut s = SiteSets::default();
        s.add_intron(201, 299); // exon ..200 | intron 201..299 | exon 300..
        assert!(s.is_donor(201));
        assert!(s.is_acceptor(299));
        assert!(!s.is_donor(299));
        assert!(!s.is_acceptor(201));
    }

    #[test]
    fn starts_and_stops_tracked() {
        let mut s = SiteSets::default();
        s.add_start(100);
        s.add_stop(400);
        assert!(s.is_start(100) && s.is_stop(400));
        assert!(!s.is_start(400));
    }
}
