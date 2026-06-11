//! GFF3 writer over the format-neutral [`OutGene`] model.

use crate::io::out_model::{OutGene, OutTranscript, cds_phases, strand_char};
use std::io::{self, Write};

/// Write genes as GFF3.
pub fn write<W: Write>(w: &mut W, genes: &[OutGene]) -> io::Result<()> {
    writeln!(w, "##gff-version 3")?;
    for gene in genes {
        let strand = strand_char(gene.strand);
        writeln!(
            w,
            "{}\tcombinr\tgene\t{}\t{}\t.\t{strand}\t.\tID={}",
            gene.contig, gene.lend, gene.rend, gene.gene_id
        )?;
        for t in &gene.transcripts {
            write_transcript(w, &gene.gene_id, t)?;
        }
    }
    Ok(())
}

fn write_transcript<W: Write>(w: &mut W, gene_id: &str, t: &OutTranscript) -> io::Result<()> {
    let strand = strand_char(t.strand);
    let (lend, rend) = (
        t.exons.iter().map(|c| c.lend).min().unwrap(),
        t.exons.iter().map(|c| c.rend).max().unwrap(),
    );
    let tid = &t.transcript_id;

    let mut attrs = format!("ID={tid};Parent={gene_id}");
    for (k, vals) in &t.attrs {
        let v = vals.iter().map(|s| escape(s)).collect::<Vec<_>>().join(",");
        attrs.push_str(&format!(";{k}={v}"));
    }
    writeln!(
        w,
        "{}\tcombinr\tmRNA\t{lend}\t{rend}\t.\t{strand}\t.\t{attrs}",
        t.contig
    )?;
    for (k, ex) in t.exons.iter().enumerate() {
        writeln!(
            w,
            "{}\tcombinr\texon\t{}\t{}\t.\t{strand}\t.\tID={tid}.exon{};Parent={tid}",
            t.contig,
            ex.lend,
            ex.rend,
            k + 1
        )?;
    }
    if !t.cds.is_empty() {
        write_utr(w, t, "five_prime_UTR", &t.five_utr)?;
        let phases = cds_phases(&t.cds, t.strand);
        for (k, seg) in t.cds.iter().enumerate() {
            let phase = phases[&(seg.lend, seg.rend)];
            writeln!(
                w,
                "{}\tcombinr\tCDS\t{}\t{}\t.\t{strand}\t{phase}\tID={tid}.cds{};Parent={tid}",
                t.contig,
                seg.lend,
                seg.rend,
                k + 1
            )?;
        }
        write_utr(w, t, "three_prime_UTR", &t.three_utr)?;
    }
    Ok(())
}

fn write_utr<W: Write>(
    w: &mut W,
    t: &OutTranscript,
    feature: &str,
    utr: &[crate::model::Coordset],
) -> io::Result<()> {
    let strand = strand_char(t.strand);
    let tid = &t.transcript_id;
    for (k, seg) in utr.iter().enumerate() {
        writeln!(
            w,
            "{}\tcombinr\t{feature}\t{}\t{}\t.\t{strand}\t.\tID={tid}.{feature}{};Parent={tid}",
            t.contig,
            seg.lend,
            seg.rend,
            k + 1
        )?;
    }
    Ok(())
}

/// Percent-encode GFF3 column-9 reserved characters in an attribute value.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            '&' => out.push_str("%26"),
            ',' => out.push_str("%2C"),
            '\t' => out.push_str("%09"),
            '\n' => out.push_str("%0A"),
            '%' => out.push_str("%25"),
            _ => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assemble::assemble_cluster;
    use crate::io::out_model::from_assemblies;
    use crate::model::{Alignment, Segment, Strand};

    fn spliced(acc: &str, contig: &str, orient: Strand, segs: &[(i64, i64)]) -> Alignment {
        let mut a = Alignment::new(
            acc,
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            orient,
        );
        a.contig = contig.to_string();
        a.spliced_orient = orient;
        a
    }

    #[test]
    fn writes_well_formed_gff3() {
        let aligns = vec![
            spliced("a", "chr1", Strand::Plus, &[(100, 200), (300, 400)]),
            spliced("b", "chr1", Strand::Plus, &[(120, 200), (300, 450)]),
        ];
        let asms = assemble_cluster(&aligns, 20).unwrap();
        let genes = from_assemblies(&asms);
        let mut buf = Vec::new();
        write(&mut buf, &genes).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.starts_with("##gff-version 3\n"));
        assert!(out.contains("chr1\tcombinr\tgene\t100\t450\t"));
        assert!(out.contains("\tmRNA\t100\t450\t"));
        assert!(out.contains("contains=a,b") || out.contains("contains=b,a"));
        assert!(out.contains("num_contained=2"));
        assert!(out.contains("\texon\t100\t200\t"));
        assert!(out.contains("\texon\t300\t450\t"));
    }
}
