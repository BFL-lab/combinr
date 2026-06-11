//! Tab-separated alternative-splicing event report.

use crate::altsplice::EventRecord;
use crate::model::Strand;
use std::io::{self, Write};

/// Write the alt-splice events as a TSV with a header row.
pub fn write_events<W: Write>(w: &mut W, events: &[EventRecord]) -> io::Result<()> {
    writeln!(
        w,
        "contig\tlocus\tevent_type\tstrand\tcoords\tisoform_a\tisoform_b\tregion\tdetail"
    )?;
    for e in events {
        let locus = e.isoform_a.split(".iso").next().unwrap_or("");
        let coords = e
            .coords
            .iter()
            .map(|c| format!("{}-{}", c.lend, c.rend))
            .collect::<Vec<_>>()
            .join(",");
        let region = e.region.map(|r| r.as_str()).unwrap_or("");
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.contig,
            locus,
            e.kind.as_str(),
            strand_char(e.strand),
            coords,
            e.isoform_a,
            e.isoform_b,
            region,
            e.detail.as_deref().unwrap_or("")
        )?;
    }
    Ok(())
}

fn strand_char(s: Strand) -> char {
    match s {
        Strand::Plus => '+',
        Strand::Minus => '-',
        Strand::Unknown => '.',
    }
}
