//! GTF writer over the format-neutral [`OutGene`] model.

use crate::io::out_model::{OutGene, OutTranscript, cds_phases, strand_char};
use std::fmt::Write as _;
use std::io::{self, Write};

/// Write genes as GTF (gene_id / transcript_id attributes, `transcript`/`exon`/
/// `CDS`/`5UTR`/`3UTR` features).
pub fn write<W: Write>(w: &mut W, genes: &[OutGene]) -> io::Result<()> {
    for gene in genes {
        let strand = strand_char(gene.strand);
        writeln!(
            w,
            "{}\tcombinr\tgene\t{}\t{}\t.\t{strand}\t.\t{}",
            gene.contig,
            gene.lend,
            gene.rend,
            gtf_attrs(&gene.gene_id, None, &[])
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
    let attrs = gtf_attrs(gene_id, Some(tid), &t.attrs);

    writeln!(
        w,
        "{}\tcombinr\ttranscript\t{lend}\t{rend}\t.\t{strand}\t.\t{attrs}",
        t.contig
    )?;
    let base = gtf_attrs(gene_id, Some(tid), &[]);
    for ex in &t.exons {
        writeln!(
            w,
            "{}\tcombinr\texon\t{}\t{}\t.\t{strand}\t.\t{base}",
            t.contig, ex.lend, ex.rend
        )?;
    }
    if !t.cds.is_empty() {
        for seg in &t.five_utr {
            writeln!(
                w,
                "{}\tcombinr\t5UTR\t{}\t{}\t.\t{strand}\t.\t{base}",
                t.contig, seg.lend, seg.rend
            )?;
        }
        let phases = cds_phases(&t.cds, t.strand, t.cds_start_phase);
        for seg in &t.cds {
            let phase = phases[&(seg.lend, seg.rend)];
            writeln!(
                w,
                "{}\tcombinr\tCDS\t{}\t{}\t.\t{strand}\t{phase}\t{base}",
                t.contig, seg.lend, seg.rend
            )?;
        }
        for seg in &t.three_utr {
            writeln!(
                w,
                "{}\tcombinr\t3UTR\t{}\t{}\t.\t{strand}\t.\t{base}",
                t.contig, seg.lend, seg.rend
            )?;
        }
    }
    Ok(())
}

/// Build a GTF column-9 attribute string. Always leads with `gene_id`, then
/// `transcript_id` (when given), then any extras.
fn gtf_attrs(
    gene_id: &str,
    transcript_id: Option<&str>,
    extra: &[(String, Vec<String>)],
) -> String {
    let mut s = format!("gene_id \"{}\";", gtf_escape(gene_id));
    if let Some(tid) = transcript_id {
        let _ = write!(s, " transcript_id \"{}\";", gtf_escape(tid));
    }
    for (k, vals) in extra {
        let v = vals
            .iter()
            .map(|x| gtf_escape(x))
            .collect::<Vec<_>>()
            .join(",");
        let _ = write!(s, " {k} \"{v}\";");
    }
    s
}

/// Escape characters that would break a GTF quoted value.
fn gtf_escape(s: &str) -> String {
    s.replace('"', "'").replace(['\t', '\n'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assemble::assemble_cluster;
    use crate::io::out_model::from_assemblies;
    use crate::model::{Alignment, Segment, Strand};

    #[test]
    fn writes_well_formed_gtf() {
        let mut a = Alignment::new(
            "a",
            vec![Segment::new(100, 200), Segment::new(300, 400)],
            Strand::Plus,
        );
        a.contig = "chr1".into();
        a.spliced_orient = Strand::Plus;
        let asms = assemble_cluster(&[a], 20).unwrap();
        let genes = from_assemblies(&asms);
        let mut buf = Vec::new();
        write(&mut buf, &genes).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("chr1\tcombinr\ttranscript\t100\t400\t"));
        assert!(out.contains("gene_id \"asm1\";"));
        assert!(out.contains("transcript_id \"asm1.t1\";"));
        assert!(out.contains("\texon\t100\t200\t"));
        assert!(out.contains("contains \"a\";"));
    }

    #[test]
    fn five_prime_partial_first_cds_row_carries_the_start_phase() {
        use crate::model::Coordset;
        let cds = vec![Coordset::new(10, 20), Coordset::new(30, 40)];
        let genes = vec![OutGene {
            gene_id: "g1".into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            lend: 10,
            rend: 40,
            attrs: vec![],
            transcripts: vec![OutTranscript {
                transcript_id: "g1.t1".into(),
                contig: "chr1".into(),
                strand: Strand::Plus,
                exons: cds.clone(),
                cds,
                five_utr: vec![],
                three_utr: vec![],
                attrs: vec![],
                cds_start_phase: 1,
            }],
        }];
        let mut buf = Vec::new();
        write(&mut buf, &genes).unwrap();
        let out = String::from_utf8(buf).unwrap();
        // 10..20 phase 1; 30..40 after 11-1 = 10 bases -> 2
        assert!(out.contains("\tCDS\t10\t20\t.\t+\t1\t"), "{out}");
        assert!(out.contains("\tCDS\t30\t40\t.\t+\t2\t"), "{out}");
    }
}
