//! Per-gene evidence report for `consensus --evidence-report` (EVidenceModeler `.evm.out`
//! style), so consensus genes can be screened by the evidence that supports them.
//!
//! Layout: two leading `##` lines describing the format, then per gene a header line
//!
//! ```text
//! # <gene_id> <contig>:<lend>-<rend> orient(<+|->) score(<%.2f>) noncoding_equivalent(<%.2f|NA>) raw_noncoding(<%.2f|NA>) S-ratio(<%.2f|inf|NA>) coding_length(<n>) low_support(<bool>) partial5(<bool>) partial3(<bool>) support(<consensus|transcript_orf>)
//! ```
//!
//! followed by one line per exon and intron in ascending genomic order, and a blank line:
//!
//! ```text
//! end5<TAB>end3<TAB><type><+|-><TAB><start_frame><TAB><end_frame><TAB>{accession;source},...
//! end5<TAB>end3<TAB>INTRON<TAB><TAB><TAB>{accession;source},...
//! ```
//!
//! `end5 > end3` on the minus strand; frames are EVM's (1-3 forward, 4-6 reverse); the
//! `source` is the evidence's GFF column 2. The gene ID is the one in the GFF3/GTF output
//! (callers pass [`crate::consensus::output::ordered_with_ids`]). Promoted transcript-ORF
//! genes print `NA` for the noncoding fields and list every contained transcript on
//! every feature (see [`crate::consensus::promote`]).

use crate::consensus::engine::{CalledGene, FeatureKind};
use crate::consensus::output::format_ratio;
use crate::model::Strand;
use std::io::{self, Write};

/// Write the report for `genes` (`(gene_id, gene)` in output order).
pub fn write_report<W: Write>(w: &mut W, genes: &[(String, &CalledGene)]) -> io::Result<()> {
    writeln!(
        w,
        "## combinr consensus evidence report (EVidenceModeler .evm.out style)"
    )?;
    writeln!(
        w,
        "## exon: end5<TAB>end3<TAB>type+strand<TAB>start_frame<TAB>end_frame<TAB>{{accession;source}},...   intron: end5<TAB>end3<TAB>INTRON<TAB><TAB><TAB>{{accession;source}},..."
    )?;
    for (gene_id, g) in genes {
        write_gene(w, gene_id, g)?;
    }
    Ok(())
}

fn write_gene<W: Write>(w: &mut W, gene_id: &str, g: &CalledGene) -> io::Result<()> {
    let (lend, rend) = g.span();
    let orient = if g.orient == Strand::Minus { '-' } else { '+' };
    let noncoding = |x: f64| {
        if g.promoted {
            "NA".to_string()
        } else {
            format!("{x:.2}")
        }
    };
    writeln!(
        w,
        "# {gene_id} {}:{lend}-{rend} orient({orient}) score({:.2}) noncoding_equivalent({}) raw_noncoding({}) S-ratio({}) coding_length({}) low_support({}) partial5({}) partial3({}) support({})",
        g.contig,
        g.score,
        noncoding(g.support.noncoding_equivalent),
        noncoding(g.support.raw_noncoding),
        format_ratio(g),
        g.support.coding_length,
        g.support.low_support,
        g.partial5,
        g.partial3,
        if g.promoted {
            "transcript_orf"
        } else {
            "consensus"
        },
    )?;
    for f in &g.features {
        let (end5, end3) = if g.orient == Strand::Minus {
            (f.coords.rend, f.coords.lend)
        } else {
            (f.coords.lend, f.coords.rend)
        };
        match f.kind {
            FeatureKind::Exon {
                exon_type,
                start_frame,
                end_frame,
            } => write!(
                w,
                "{end5}\t{end3}\t{}{orient}\t{start_frame}\t{end_frame}\t",
                exon_type.as_str()
            )?,
            FeatureKind::Intron => write!(w, "{end5}\t{end3}\tINTRON\t\t\t")?,
        }
        let mut seen: Vec<&(String, String)> = Vec::with_capacity(f.evidence.len());
        for ev in &f.evidence {
            if !seen.contains(&ev) {
                seen.push(ev);
            }
        }
        let tokens: Vec<String> = seen.iter().map(|(a, s)| format!("{{{a};{s}}}")).collect();
        writeln!(w, "{}", tokens.join(","))?;
    }
    writeln!(w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::engine::FeatureSupport;
    use crate::consensus::exon::ExonType;
    use crate::consensus::filter::SupportFlags;
    use crate::model::Coordset;

    fn ev(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|&(a, s)| (a.to_string(), s.to_string()))
            .collect()
    }

    fn exon(l: i64, r: i64, t: ExonType, sf: u8, ef: u8, e: &[(&str, &str)]) -> FeatureSupport {
        FeatureSupport {
            coords: Coordset::new(l, r),
            kind: FeatureKind::Exon {
                exon_type: t,
                start_frame: sf,
                end_frame: ef,
            },
            evidence: ev(e),
        }
    }

    fn intron(l: i64, r: i64, e: &[(&str, &str)]) -> FeatureSupport {
        FeatureSupport {
            coords: Coordset::new(l, r),
            kind: FeatureKind::Intron,
            evidence: ev(e),
        }
    }

    fn gene(orient: Strand, promoted: bool, features: Vec<FeatureSupport>) -> CalledGene {
        let exons: Vec<Coordset> = features
            .iter()
            .filter(|f| matches!(f.kind, FeatureKind::Exon { .. }))
            .map(|f| f.coords)
            .collect();
        CalledGene {
            contig: "chr1".into(),
            orient,
            exons: exons.clone(),
            cds: exons,
            five_utr: vec![],
            three_utr: vec![],
            partial5: false,
            partial3: true,
            cds_start_phase: 0,
            score: 1234.567,
            support: SupportFlags {
                raw_noncoding: 0.0,
                noncoding_equivalent: if promoted { 0.0 } else { 0.1235 },
                score_ratio: if promoted { f64::INFINITY } else { 9996.49 },
                coding_length: 60,
                low_support: false,
            },
            promoted,
            features,
        }
    }

    fn render(genes: &[(String, &CalledGene)]) -> Vec<String> {
        let mut buf = Vec::new();
        write_report(&mut buf, genes).unwrap();
        String::from_utf8(buf)
            .unwrap()
            .lines()
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn plus_and_minus_genes_render_evm_style() {
        let plus = gene(
            Strand::Plus,
            false,
            vec![
                exon(
                    10,
                    39,
                    ExonType::Initial,
                    1,
                    3,
                    &[("p1", "fgenesh"), ("a1", "gmap")],
                ),
                intron(40, 99, &[("a1", "gmap"), ("a1", "gmap"), ("p1", "fgenesh")]),
                exon(100, 129, ExonType::Terminal, 1, 3, &[("p1", "fgenesh")]),
            ],
        );
        let minus = gene(
            Strand::Minus,
            false,
            vec![
                exon(500, 529, ExonType::Terminal, 4, 6, &[("q", "genewise")]),
                intron(530, 599, &[("q", "genewise")]),
                exon(600, 629, ExonType::Initial, 5, 4, &[("q", "genewise")]),
            ],
        );
        let lines = render(&[
            ("consensus.chr1.g1".into(), &plus),
            ("consensus.chr1.g2".into(), &minus),
        ]);
        let want = [
            "# consensus.chr1.g1 chr1:10-129 orient(+) score(1234.57) noncoding_equivalent(0.12) raw_noncoding(0.00) S-ratio(9996.49) coding_length(60) low_support(false) partial5(false) partial3(true) support(consensus)",
            "10\t39\tinitial+\t1\t3\t{p1;fgenesh},{a1;gmap}",
            "40\t99\tINTRON\t\t\t{a1;gmap},{p1;fgenesh}",
            "100\t129\tterminal+\t1\t3\t{p1;fgenesh}",
            "",
            "# consensus.chr1.g2 chr1:500-629 orient(-) score(1234.57) noncoding_equivalent(0.12) raw_noncoding(0.00) S-ratio(9996.49) coding_length(60) low_support(false) partial5(false) partial3(true) support(consensus)",
            "529\t500\tterminal-\t4\t6\t{q;genewise}",
            "599\t530\tINTRON\t\t\t{q;genewise}",
            "629\t600\tinitial-\t5\t4\t{q;genewise}",
            "",
        ];
        assert!(lines[0].starts_with("## ") && lines[1].starts_with("## "));
        assert_eq!(&lines[2..], &want[..]);
    }

    #[test]
    fn promoted_gene_prints_na_noncoding_fields() {
        let g = gene(
            Strand::Plus,
            true,
            vec![exon(
                10,
                69,
                ExonType::Single,
                1,
                3,
                &[("t1", "PASA"), ("t2", "transcript")],
            )],
        );
        let lines = render(&[("consensus.chr1.g1".into(), &g)]);
        assert_eq!(
            lines[2],
            "# consensus.chr1.g1 chr1:10-69 orient(+) score(1234.57) noncoding_equivalent(NA) raw_noncoding(NA) S-ratio(NA) coding_length(60) low_support(false) partial5(false) partial3(true) support(transcript_orf)"
        );
        assert_eq!(lines[3], "10\t69\tsingle+\t1\t3\t{t1;PASA},{t2;transcript}");
        assert_eq!(lines[4], "");
        assert_eq!(lines.len(), 5);
    }
}
