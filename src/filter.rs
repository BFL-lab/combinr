//! Optional pre-assembly quality filters. All thresholds are **off by default**
//! (combinr trusts its input); enable them to police raw aligner output.
//!
//! `min_percent_aligned` from PASA is intentionally absent: it needs the full
//! transcript length, which GTF/GFF3 alignments do not carry.

use crate::model::Alignment;

/// Quality thresholds. `None` means "do not filter on this".
#[derive(Debug, Clone, Default)]
pub struct Filters {
    /// Minimum length-weighted average percent identity (when present on input).
    pub min_avg_per_id: Option<f64>,
    /// Minimum intron length (bp); shorter introns invalidate the alignment.
    pub min_intron: Option<i64>,
    /// Maximum intron length (bp); longer introns invalidate the alignment.
    pub max_intron: Option<i64>,
}

impl Filters {
    pub fn none() -> Self {
        Filters::default()
    }

    pub fn is_noop(&self) -> bool {
        self.min_avg_per_id.is_none() && self.min_intron.is_none() && self.max_intron.is_none()
    }

    /// `true` if `a` passes all enabled thresholds.
    pub fn keeps(&self, a: &Alignment) -> bool {
        if let Some(min) = self.min_avg_per_id
            && let Some(pid) = a.per_id
            && pid < min
        {
            return false;
        }
        if self.min_intron.is_some() || self.max_intron.is_some() {
            for w in a.segments.windows(2) {
                let intron_len = (w[1].coords.lend - 1) - (w[0].coords.rend + 1) + 1;
                if let Some(mn) = self.min_intron
                    && intron_len < mn
                {
                    return false;
                }
                if let Some(mx) = self.max_intron
                    && intron_len > mx
                {
                    return false;
                }
            }
        }
        true
    }
}

/// Drop alignments that fail any enabled threshold. Returns the kept alignments
/// and the number dropped (for reporting).
pub fn apply(alignments: Vec<Alignment>, filters: &Filters) -> (Vec<Alignment>, usize) {
    if filters.is_noop() {
        return (alignments, 0);
    }
    let total = alignments.len();
    let kept: Vec<Alignment> = alignments
        .into_iter()
        .filter(|a| filters.keeps(a))
        .collect();
    let dropped = total - kept.len();
    (kept, dropped)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Segment, Strand};

    fn al(per_id: Option<f64>, segs: &[(i64, i64)]) -> Alignment {
        let mut a = Alignment::new(
            "x",
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            Strand::Plus,
        );
        a.per_id = per_id;
        a
    }

    #[test]
    fn noop_keeps_everything() {
        let f = Filters::none();
        assert!(f.is_noop());
        let (kept, dropped) = apply(vec![al(Some(50.0), &[(1, 10)])], &f);
        assert_eq!(kept.len(), 1);
        assert_eq!(dropped, 0);
    }

    #[test]
    fn per_id_threshold() {
        let f = Filters {
            min_avg_per_id: Some(95.0),
            ..Default::default()
        };
        assert!(f.keeps(&al(Some(99.0), &[(1, 10)])));
        assert!(!f.keeps(&al(Some(90.0), &[(1, 10)])));
        // missing per_id is not filtered out (cannot evaluate)
        assert!(f.keeps(&al(None, &[(1, 10)])));
    }

    #[test]
    fn intron_length_bounds() {
        // exons 1-10 and 31-40 → intron 11..30 = length 20
        let a = al(None, &[(1, 10), (31, 40)]);
        let short = Filters {
            min_intron: Some(21),
            ..Default::default()
        };
        assert!(!short.keeps(&a), "intron of 20 < min 21 is rejected");
        let long = Filters {
            max_intron: Some(19),
            ..Default::default()
        };
        assert!(!long.keeps(&a), "intron of 20 > max 19 is rejected");
        let ok = Filters {
            min_intron: Some(10),
            max_intron: Some(100),
            ..Default::default()
        };
        assert!(ok.keeps(&a));
    }
}
