//! The five alternative-splicing classifiers, ported from
//! `Alternative_splice_comparer.pm`. Every overlap test is **strict**
//! (`overlaps_strict`), matching the Perl `i_lend < j_rend && i_rend > j_lend`.

use super::{EventKind, EventRecord, Isoform};
use crate::model::{Coordset, Strand};

/// Compare two isoforms and return all alternative-splicing events between them,
/// running each classifier in the same directions as the Perl driver.
pub fn compare_isoforms(a: &Isoform, b: &Isoform, fuzz: i64) -> Vec<EventRecord> {
    let mut out = Vec::new();

    // retained introns: both directions
    out.extend(find_retained_introns(a, b));
    out.extend(find_retained_introns(b, a));

    // conventional alt donor/acceptor: once (internally symmetric)
    out.extend(find_conventional(a, b));

    // transcription start / polyA within intron: both directions
    out.extend(find_starts_and_ends_within_introns(a, b, fuzz));
    out.extend(find_starts_and_ends_within_introns(b, a, fuzz));

    // exon skipping: both directions
    out.extend(find_exon_skips(a, b));
    out.extend(find_exon_skips(b, a));

    // alternate terminal exons: both directions
    out.extend(find_alternate_exons(a, b));
    out.extend(find_alternate_exons(b, a));

    out
}

fn record(
    kind: EventKind,
    a: &Isoform,
    b: &Isoform,
    coords: Vec<Coordset>,
    detail: Option<String>,
) -> EventRecord {
    EventRecord {
        kind,
        contig: a.contig.clone(),
        strand: a.strand,
        coords,
        isoform_a: a.id.clone(),
        isoform_b: b.id.clone(),
        detail,
        region: None,
    }
}

/// Retained (unspliced) intron: an intron of `b` lies entirely within an exon of
/// `a`. (`find_unspliced_introns`.)
pub fn find_retained_introns(a: &Isoform, b: &Isoform) -> Vec<EventRecord> {
    let mut out = Vec::new();
    for intron in b.introns() {
        for exon in &a.exons {
            if intron.lend > exon.lend && intron.rend < exon.rend {
                out.push(record(EventKind::RetainedIntron, a, b, vec![intron], None));
            }
        }
    }
    out
}

/// Conventional alt-splice: differing donor/acceptor sites at 1-to-1 mapped
/// internal exons. (`find_conventional_alt_splice_isoforms`.) Runs once.
pub fn find_conventional(a: &Isoform, b: &Isoform) -> Vec<EventRecord> {
    let ax = &a.exons; // lend-sorted
    let bx = &b.exons;
    let last_a = ax.len().saturating_sub(1);
    let last_b = bx.len().saturating_sub(1);

    // all-vs-all strict overlap → match indices
    let mut a_match: Vec<Vec<usize>> = vec![Vec::new(); ax.len()];
    let mut b_match: Vec<Vec<usize>> = vec![Vec::new(); bx.len()];
    for (i, ae) in ax.iter().enumerate() {
        for (j, be) in bx.iter().enumerate() {
            if ae.overlaps_strict(be) {
                a_match[i].push(j);
                b_match[j].push(i);
            }
        }
    }

    let plus = a.strand != Strand::Minus;
    let mut out = Vec::new();
    for i in 0..ax.len() {
        if a_match[i].len() != 1 {
            continue;
        }
        let j = a_match[i][0];
        if b_match[j].len() != 1 || b_match[j][0] != i {
            continue; // require mutual 1-to-1
        }
        // left boundary difference at an internal junction
        if ax[i].lend != bx[j].lend && i != 0 && j != 0 {
            let kind = if plus {
                EventKind::AltAcceptor
            } else {
                EventKind::AltDonor
            };
            out.push(record(
                kind,
                a,
                b,
                vec![Coordset::new(ax[i].lend, bx[j].lend)],
                None,
            ));
        }
        // right boundary difference at an internal junction
        if ax[i].rend != bx[j].rend && i != last_a && j != last_b {
            let kind = if plus {
                EventKind::AltDonor
            } else {
                EventKind::AltAcceptor
            };
            out.push(record(
                kind,
                a,
                b,
                vec![Coordset::new(ax[i].rend, bx[j].rend)],
                None,
            ));
        }
    }
    out
}

/// Transcription start / polyadenylation within an intron of the other isoform.
/// (`find_starts_and_ends_within_introns`.)
// The `abs(...) + 1 <= fuzz` distance checks mirror the Perl verbatim; clippy
// would rewrite them to `< fuzz` (mathematically identical for integers) but the
// `+ 1` keeps the port readable against the source.
#[allow(clippy::int_plus_one)]
pub fn find_starts_and_ends_within_introns(
    a: &Isoform,
    b: &Isoform,
    fuzz: i64,
) -> Vec<EventRecord> {
    let mut out = Vec::new();
    let plus = a.strand != Strand::Minus;
    let b_introns = b.introns();

    // start: a's 5' terminus inside a b-intron
    let e5 = a.tss_coord();
    let tss_exon = a.first_transcript_exon();
    for intron in &b_introns {
        if e5 >= intron.lend && e5 <= intron.rend {
            let is_fuzz = if plus {
                (e5 - intron.rend).abs() + 1 <= fuzz
            } else {
                (e5 - intron.lend).abs() + 1 <= fuzz
            };
            if is_fuzz {
                continue;
            }
            if b.exons.iter().any(|ex| tss_exon.overlaps_strict(ex)) {
                out.push(record(
                    EventKind::StartWithinIntron,
                    a,
                    b,
                    vec![Coordset::new(e5, e5)],
                    None,
                ));
            }
            break;
        }
    }

    // end: a's 3' terminus inside a b-intron
    let e3 = a.tes_coord();
    let tes_exon = a.last_transcript_exon();
    for intron in &b_introns {
        if e3 >= intron.lend && e3 <= intron.rend {
            let is_fuzz = if plus {
                (e3 - intron.lend).abs() + 1 <= fuzz
            } else {
                (e3 - intron.rend).abs() + 1 <= fuzz
            };
            if is_fuzz {
                continue;
            }
            if b.exons.iter().any(|ex| tes_exon.overlaps_strict(ex)) {
                out.push(record(
                    EventKind::EndWithinIntron,
                    a,
                    b,
                    vec![Coordset::new(e3, e3)],
                    None,
                ));
            }
            break;
        }
    }

    out
}

/// Exon skipping: an internal exon of `a` lies within an intron of `b`, with
/// flanking exons anchorable on both sides. Adjacent skipped exons are grouped
/// into one event. (`find_exon_skipping_events`.)
pub fn find_exon_skips(a: &Isoform, b: &Isoform) -> Vec<EventRecord> {
    let ax = &a.exons; // lend-sorted
    let b_introns = b.introns();

    // candidate skipped exons: strictly inside some b-intron
    let mut skipped: Vec<(usize, Coordset)> = Vec::new();
    for (idx, exon) in ax.iter().enumerate() {
        let inside = b_introns
            .iter()
            .any(|intron| exon.lend > intron.lend && exon.rend < intron.rend);
        if !inside {
            continue;
        }
        // left-anchorable: some a-exon fully left of this one overlaps a b-exon
        // also fully left of this one.
        let left_ok = ax.iter().any(|ae| {
            ae.rend < exon.lend
                && b.exons
                    .iter()
                    .any(|be| be.rend < exon.lend && ae.overlaps_strict(be))
        });
        // right-anchorable: mirror.
        let right_ok = ax.iter().any(|ae| {
            ae.lend > exon.rend
                && b.exons
                    .iter()
                    .any(|be| be.lend > exon.rend && ae.overlaps_strict(be))
        });
        if left_ok && right_ok {
            skipped.push((idx, *exon));
        }
    }

    // group exons adjacent in a's exon order into single events.
    let mut out = Vec::new();
    let mut group: Vec<Coordset> = Vec::new();
    let mut prev_idx: Option<usize> = None;
    for (idx, exon) in skipped {
        match prev_idx {
            Some(p) if idx == p + 1 => group.push(exon),
            _ => {
                if !group.is_empty() {
                    out.push(record(
                        EventKind::ExonSkip,
                        a,
                        b,
                        std::mem::take(&mut group),
                        None,
                    ));
                }
                group.push(exon);
            }
        }
        prev_idx = Some(idx);
    }
    if !group.is_empty() {
        out.push(record(EventKind::ExonSkip, a, b, group, None));
    }
    out
}

/// Alternate terminal exons: non-overlapping terminal exons of `a` preceding the
/// first (front) / following the last (back) exon that overlaps `b`.
/// (`find_alternate_exons`.) "front"/"back" are genomic sides.
pub fn find_alternate_exons(a: &Isoform, b: &Isoform) -> Vec<EventRecord> {
    let ax = &a.exons; // lend-sorted
    let bx = &b.exons;
    let last_a = ax.len().saturating_sub(1);
    let last_b = bx.len().saturating_sub(1);
    let mut out = Vec::new();

    // front: first a-exon overlapping any b-exon
    let mut front: Option<(usize, usize)> = None;
    'f: for (i, ae) in ax.iter().enumerate() {
        for (j, be) in bx.iter().enumerate() {
            if ae.overlaps_strict(be) {
                front = Some((i, j));
                break 'f;
            }
        }
    }
    if let Some((i, j)) = front
        && i != 0
        && j != 0
    {
        let alt: Vec<Coordset> = ax[..i].to_vec();
        if let Some(region) = region_of(&alt) {
            out.push(record(
                EventKind::AlternateExon,
                a,
                b,
                vec![region],
                Some(format!("side=front;num_exons={}", alt.len())),
            ));
        }
    }

    // back: last a-exon overlapping any b-exon
    let mut back: Option<(usize, usize)> = None;
    'b: for i in (0..ax.len()).rev() {
        for j in (0..bx.len()).rev() {
            if ax[i].overlaps_strict(&bx[j]) {
                back = Some((i, j));
                break 'b;
            }
        }
    }
    if let Some((i, j)) = back
        && i != last_a
        && j != last_b
    {
        let alt: Vec<Coordset> = ax[(i + 1)..].to_vec();
        if let Some(region) = region_of(&alt) {
            out.push(record(
                EventKind::AlternateExon,
                a,
                b,
                vec![region],
                Some(format!("side=back;num_exons={}", alt.len())),
            ));
        }
    }

    out
}

/// Bounding region (min lend .. max rend) of a set of exons.
fn region_of(exons: &[Coordset]) -> Option<Coordset> {
    let lend = exons.iter().map(|e| e.lend).min()?;
    let rend = exons.iter().map(|e| e.rend).max()?;
    Some(Coordset { lend, rend })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn iso(id: &str, strand: Strand, exons: &[(i64, i64)]) -> Isoform {
        Isoform {
            id: id.to_string(),
            contig: "chr1".to_string(),
            strand,
            exons: exons.iter().map(|&(l, r)| Coordset::new(l, r)).collect(),
            contained_accs: vec![id.to_string()],
            source_set: BTreeSet::new(),
            is_fl: false,
        }
    }

    #[test]
    fn retained_intron_detected_one_direction() {
        // a has a single long exon covering b's intron.
        let a = iso("a", Strand::Plus, &[(100, 600)]);
        let b = iso("b", Strand::Plus, &[(100, 200), (400, 600)]);
        let in_a = find_retained_introns(&a, &b);
        assert_eq!(in_a.len(), 1);
        assert_eq!(in_a[0].kind, EventKind::RetainedIntron);
        assert_eq!(in_a[0].coords[0], Coordset::new(201, 399));
        // reverse direction: b has no intron inside an a-exon-of-b → none.
        assert!(find_retained_introns(&b, &a).is_empty());
    }

    #[test]
    fn alt_acceptor_and_donor_swap_by_strand() {
        // shared first/last exons; middle exon's left boundary differs.
        let a = iso("a", Strand::Plus, &[(100, 200), (300, 400), (500, 600)]);
        let b = iso("b", Strand::Plus, &[(100, 200), (320, 400), (500, 600)]);
        let ev = find_conventional(&a, &b);
        assert_eq!(ev.len(), 1);
        // left boundary differs on '+' → alt acceptor
        assert_eq!(ev[0].kind, EventKind::AltAcceptor);

        // same structures but minus strand → the left-boundary diff is a donor
        let am = iso("a", Strand::Minus, &[(100, 200), (300, 400), (500, 600)]);
        let bm = iso("b", Strand::Minus, &[(100, 200), (320, 400), (500, 600)]);
        let evm = find_conventional(&am, &bm);
        assert_eq!(evm.len(), 1);
        assert_eq!(evm[0].kind, EventKind::AltDonor);
    }

    #[test]
    fn exon_skip_with_anchors_and_grouping() {
        // b skips a's middle exon; flanking exons overlap b's flanks.
        let a = iso("a", Strand::Plus, &[(100, 200), (300, 400), (500, 600)]);
        let b = iso("b", Strand::Plus, &[(100, 210), (490, 600)]);
        let ev = find_exon_skips(&a, &b);
        assert_eq!(ev.len(), 1);
        assert_eq!(ev[0].kind, EventKind::ExonSkip);
        assert_eq!(ev[0].coords, vec![Coordset::new(300, 400)]);
    }

    #[test]
    fn start_within_intron_respects_fuzz() {
        // a's first exon starts at 350 (inside b's intron 201..449) and extends
        // into b's exon (450..600), so the start sits within the intron and the
        // exon still overlaps another b-exon, as the classifier requires.
        let a = iso("a", Strand::Plus, &[(350, 500), (700, 800)]);
        let b = iso("b", Strand::Plus, &[(100, 200), (450, 600)]);
        let ev = find_starts_and_ends_within_introns(&a, &b, 20);
        assert!(ev.iter().any(|e| e.kind == EventKind::StartWithinIntron));

        // If the start sits within `fuzz` of the intron's 3' boundary it is
        // treated as boundary noise and not reported.
        let near = iso("near", Strand::Plus, &[(440, 500), (700, 800)]);
        let ev2 = find_starts_and_ends_within_introns(&near, &b, 20);
        assert!(!ev2.iter().any(|e| e.kind == EventKind::StartWithinIntron));
    }

    #[test]
    fn alternate_front_exon_detected() {
        // The front rule needs the first overlapping a-exon index != 0 AND its
        // b-partner index != 0, so both isoforms must have a (different) leading
        // exon before the first shared one.
        let a = iso("a", Strand::Plus, &[(10, 50), (100, 200), (300, 400)]);
        let b = iso("b", Strand::Plus, &[(60, 90), (100, 200), (300, 400)]);
        let ev = find_alternate_exons(&a, &b);
        assert!(ev.iter().any(|e| e.kind == EventKind::AlternateExon
            && e.detail.as_deref() == Some("side=front;num_exons=1")));
    }
}
