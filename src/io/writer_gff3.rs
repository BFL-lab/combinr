//! GFF3 writer over the format-neutral [`OutGene`] model.

use crate::io::out_model::{OutGene, OutTranscript, cds_phases, strand_char};
use std::fmt::Write as _;
use std::io::{self, Write};

/// Write genes as GFF3.
pub fn write<W: Write>(w: &mut W, genes: &[OutGene]) -> io::Result<()> {
    writeln!(w, "##gff-version 3")?;
    for gene in genes {
        let strand = strand_char(gene.strand);
        let mut attrs = format!("ID={}", gene.gene_id);
        push_attrs(&mut attrs, &gene.attrs);
        writeln!(
            w,
            "{}\tcombinr\tgene\t{}\t{}\t.\t{strand}\t.\t{attrs}",
            gene.contig, gene.lend, gene.rend
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
    push_attrs(&mut attrs, &t.attrs);
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
        let phases = cds_phases(&t.cds, t.strand, t.cds_start_phase);
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

/// Append `;key=v1,v2` for each attribute, values escaped (the inverse of [`unescape`]).
fn push_attrs(out: &mut String, attrs: &[(String, Vec<String>)]) {
    for (k, vals) in attrs {
        let v = vals.iter().map(|s| escape(s)).collect::<Vec<_>>().join(",");
        let _ = write!(out, ";{k}={v}");
    }
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

/// Decode a percent-encoded GFF3 column-9 value: the exact inverse of [`escape`].
/// Every `%XX` (two hex digits) decodes to its byte; a `%` not followed by two hex
/// digits is kept literally.
pub(crate) fn unescape(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let b = s.as_bytes();
    let hex = |c: u8| (c as char).to_digit(16);
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%'
            && i + 2 < b.len()
            && let (Some(hi), Some(lo)) = (hex(b[i + 1]), hex(b[i + 2]))
        {
            out.push((hi * 16 + lo) as u8);
            i += 3;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
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

    #[test]
    fn unescape_inverts_escape() {
        for v in [
            "plain",
            "a;b=c&d,e\tf\ng%h",
            "%3B literal",
            "100%",
            "%zz",
            "ünïcode;",
            "",
        ] {
            assert_eq!(unescape(&escape(v)), v, "round-trip of {v:?}");
        }
        assert_eq!(unescape("a%3Bb%2Cc%25"), "a;b,c%");
        assert_eq!(unescape("100%"), "100%");
    }

    #[test]
    fn gene_attrs_follow_id() {
        let genes = vec![OutGene {
            gene_id: "g1".into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            lend: 1,
            rend: 10,
            attrs: vec![("Name".into(), vec!["x;y".into(), "z".into()])],
            transcripts: vec![],
        }];
        let mut buf = Vec::new();
        write(&mut buf, &genes).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("\tgene\t1\t10\t.\t+\t.\tID=g1;Name=x%3By,z\n"));
    }

    #[test]
    fn five_prime_partial_first_cds_row_carries_the_start_phase() {
        use crate::model::Coordset;
        let cds = vec![
            Coordset::new(10, 20),
            Coordset::new(30, 40),
            Coordset::new(50, 58),
        ];
        let gene = |strand: Strand| OutGene {
            gene_id: "g1".into(),
            contig: "chr1".into(),
            strand,
            lend: 10,
            rend: 58,
            attrs: vec![],
            transcripts: vec![OutTranscript {
                transcript_id: "g1.t1".into(),
                contig: "chr1".into(),
                strand,
                exons: cds.clone(),
                cds: cds.clone(),
                five_utr: vec![],
                three_utr: vec![],
                attrs: vec![],
                cds_start_phase: 2,
            }],
        };
        let render = |strand| {
            let mut buf = Vec::new();
            write(&mut buf, &[gene(strand)]).unwrap();
            String::from_utf8(buf).unwrap()
        };
        // plus, 5'->3': 10..20 phase 2; 30..40 after 11-2 = 9 bases -> 0; 50..58 after
        // 20 -> 1
        let plus = render(Strand::Plus);
        assert!(plus.contains("\tCDS\t10\t20\t.\t+\t2\t"), "{plus}");
        assert!(plus.contains("\tCDS\t30\t40\t.\t+\t0\t"), "{plus}");
        assert!(plus.contains("\tCDS\t50\t58\t.\t+\t1\t"), "{plus}");
        assert!(!plus.contains("five_prime_UTR"));
        // minus, 5'->3': 50..58 phase 2; 30..40 after 9-2 = 7 -> 2; 10..20 after 18 -> 0
        let minus = render(Strand::Minus);
        assert!(minus.contains("\tCDS\t50\t58\t.\t-\t2\t"), "{minus}");
        assert!(minus.contains("\tCDS\t30\t40\t.\t-\t2\t"), "{minus}");
        assert!(minus.contains("\tCDS\t10\t20\t.\t-\t0\t"), "{minus}");
    }
}
