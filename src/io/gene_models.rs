//! Gene-model GFF3 reader: a full `gene` → `mRNA`/`transcript` → `exon`/`CDS`/
//! `five_prime_UTR`/`three_prime_UTR` hierarchy (linked by `Parent=`) into the
//! format-neutral [`OutGene`] model, so an existing gene set can be re-emitted as-is
//! (`assemble --models`).
//!
//! - Genes come out in first-appearance order (of the gene row or of any row belonging
//!   to it), transcripts in first-appearance order under their gene.
//! - An mRNA whose `Parent` has no gene row gets a synthesized gene with that id; an mRNA
//!   with no `Parent` gets a synthesized gene `{mRNA id}.gene`. A child row whose
//!   `Parent` is no mRNA row gets a synthesized mRNA: with id `{gene}.mRNA` under that
//!   gene when the parent is a gene row, else with the parent's id under a synthesized
//!   gene `{parent}.gene`. (Synthesized ids never reuse an input id: GFF3 ids are unique.)
//! - exons/CDS/UTRs are lend-sorted; an mRNA with no exon rows gets exons = the coalesced
//!   union of its CDS and UTRs (or its own span when it has none of these either).
//!   `cds_start_phase` is the column-8 phase of the 5'-most CDS row (`.` → 0).
//! - Attributes are read in file order and percent-decoded (the inverse of the GFF3
//!   writer's escaping), each value split on `,`. Every gene attribute but `ID`, and every
//!   mRNA attribute but `ID`/`Parent`, is kept. `ID`/`Parent` values are kept verbatim.
//! - Rows of any other feature type (and children parented only by such features, e.g. a
//!   tRNA's exons) are ignored.

use crate::error::{CombinrError, Result};
use crate::io::gff;
use crate::io::out_model::{OutGene, OutTranscript};
use crate::io::writer_gff3::unescape;
use crate::model::{Coordset, Strand, merge_coords};
use std::collections::{HashMap, HashSet};
use std::path::Path;

type Attrs = Vec<(String, Vec<String>)>;

/// Column 9 in file order: `(key, raw value)` pairs, keys and values trimmed.
fn raw_attrs(col9: &str) -> Vec<(&str, &str)> {
    col9.split(';')
        .filter_map(|part| {
            let (k, v) = part.trim().split_once('=')?;
            Some((k.trim(), v.trim()))
        })
        .collect()
}

/// The kept (non-`skip`) attributes, values split on `,` then percent-decoded.
fn kept_attrs(raw: &[(&str, &str)], skip: &[&str]) -> Attrs {
    raw.iter()
        .filter(|(k, _)| !skip.contains(k))
        .map(|&(k, v)| (k.to_string(), v.split(',').map(unescape).collect()))
        .collect()
}

fn attr<'a>(raw: &[(&'a str, &'a str)], key: &str) -> Option<&'a str> {
    raw.iter().find(|(k, _)| *k == key).map(|&(_, v)| v)
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Gene,
    Mrna,
    Exon,
    Cds,
    FiveUtr,
    ThreeUtr,
}

struct Row<'a> {
    line: usize,
    kind: Kind,
    contig: &'a str,
    coords: Coordset,
    strand: Strand,
    phase: u8,
    attrs: Vec<(&'a str, &'a str)>,
}

struct TxAcc {
    id: String,
    contig: String,
    strand: Strand,
    /// The mRNA row's own span, when there is a row.
    span: Option<Coordset>,
    attrs: Attrs,
    exons: Vec<Coordset>,
    cds: Vec<(Coordset, u8)>,
    five: Vec<Coordset>,
    three: Vec<Coordset>,
}

struct GeneAcc {
    id: String,
    contig: String,
    strand: Strand,
    span: Option<Coordset>,
    attrs: Attrs,
    tx: Vec<usize>,
}

/// Load a gene-model GFF3 file (see the module doc).
pub fn load_gene_models(path: &Path) -> Result<Vec<OutGene>> {
    let file = path.display().to_string();
    let text = std::fs::read_to_string(path).map_err(|e| CombinrError::Parse {
        file: file.clone(),
        line: 0,
        msg: format!("cannot read gene-model GFF3: {e}"),
    })?;
    parse_gene_models(&text, &file)
}

/// Parse gene-model GFF3 text (see the module doc). `file` labels parse errors.
pub fn parse_gene_models(text: &str, file: &str) -> Result<Vec<OutGene>> {
    let bad = |line: usize, msg: String| CombinrError::Parse {
        file: file.to_string(),
        line,
        msg,
    };

    // Pass 1: tokenize the rows we read; note the ids of gene / mRNA / other features.
    let mut rows: Vec<Row> = Vec::new();
    let mut gene_ids: HashSet<&str> = HashSet::new();
    let mut mrna_parent: HashMap<&str, Option<&str>> = HashMap::new();
    let mut other_ids: HashSet<&str> = HashSet::new();
    for (i, raw) in text.lines().enumerate() {
        let line = i + 1;
        let Some(rec) = gff::record(raw, file, line)? else {
            continue;
        };
        let cols = gff::columns(raw).expect("record() accepted the line");
        let attrs = raw_attrs(rec.attrs);
        let kind = match rec.ftype {
            "gene" => Kind::Gene,
            "mRNA" | "transcript" => Kind::Mrna,
            "exon" => Kind::Exon,
            "CDS" => Kind::Cds,
            "five_prime_UTR" => Kind::FiveUtr,
            "three_prime_UTR" => Kind::ThreeUtr,
            _ => {
                if let Some(id) = attr(&attrs, "ID") {
                    other_ids.insert(id);
                }
                continue;
            }
        };
        if !matches!(cols[6], "+" | "-" | "." | "?") {
            return Err(bad(line, format!("bad strand {:?}", cols[6])));
        }
        let phase = match cols[7] {
            "." | "0" => 0,
            "1" => 1,
            "2" => 2,
            p => return Err(bad(line, format!("bad phase {p:?}"))),
        };
        if rec.lend > rec.rend {
            return Err(bad(line, "start coordinate > end coordinate".into()));
        }
        let first_parent = attr(&attrs, "Parent").and_then(|p| p.split(',').next());
        match kind {
            Kind::Gene | Kind::Mrna => {
                let id = attr(&attrs, "ID")
                    .ok_or_else(|| bad(line, format!("{} row without ID", rec.ftype)))?;
                if kind == Kind::Gene {
                    gene_ids.insert(id);
                } else {
                    mrna_parent.insert(id, first_parent);
                }
            }
            _ if first_parent.is_none() => {
                return Err(bad(line, format!("{} row without Parent", rec.ftype)));
            }
            _ => {}
        }
        rows.push(Row {
            line,
            kind,
            contig: rec.contig,
            coords: Coordset::new(rec.lend, rec.rend),
            strand: rec.strand,
            phase,
            attrs,
        });
    }

    // Pass 2 (file order): attach each row to its gene / transcript.
    let mut genes: Vec<GeneAcc> = Vec::new();
    let mut gene_idx: HashMap<String, usize> = HashMap::new();
    let mut txs: Vec<TxAcc> = Vec::new();
    let mut tx_idx: HashMap<String, usize> = HashMap::new();

    fn gene_of(
        genes: &mut Vec<GeneAcc>,
        gene_idx: &mut HashMap<String, usize>,
        id: &str,
        row: &Row,
    ) -> usize {
        *gene_idx.entry(id.to_string()).or_insert_with(|| {
            genes.push(GeneAcc {
                id: id.to_string(),
                contig: row.contig.to_string(),
                strand: row.strand,
                span: None,
                attrs: Vec::new(),
                tx: Vec::new(),
            });
            genes.len() - 1
        })
    }
    fn tx_of(
        genes: &mut [GeneAcc],
        txs: &mut Vec<TxAcc>,
        tx_idx: &mut HashMap<String, usize>,
        gi: usize,
        id: &str,
        row: &Row,
    ) -> usize {
        *tx_idx.entry(id.to_string()).or_insert_with(|| {
            txs.push(TxAcc {
                id: id.to_string(),
                contig: row.contig.to_string(),
                strand: row.strand,
                span: None,
                attrs: Vec::new(),
                exons: Vec::new(),
                cds: Vec::new(),
                five: Vec::new(),
                three: Vec::new(),
            });
            genes[gi].tx.push(txs.len() - 1);
            txs.len() - 1
        })
    }

    for row in &rows {
        match row.kind {
            Kind::Gene => {
                let id = attr(&row.attrs, "ID").unwrap();
                let gi = gene_of(&mut genes, &mut gene_idx, id, row);
                let g = &mut genes[gi];
                if g.span.is_some() {
                    return Err(bad(row.line, format!("duplicate gene ID {id:?}")));
                }
                (g.contig, g.strand) = (row.contig.to_string(), row.strand);
                g.span = Some(row.coords);
                g.attrs = kept_attrs(&row.attrs, &["ID"]);
            }
            Kind::Mrna => {
                let id = attr(&row.attrs, "ID").unwrap();
                let gene_id = match mrna_parent[id] {
                    Some(p) => p.to_string(),
                    None => format!("{id}.gene"),
                };
                let gi = gene_of(&mut genes, &mut gene_idx, &gene_id, row);
                let ti = tx_of(&mut genes, &mut txs, &mut tx_idx, gi, id, row);
                let t = &mut txs[ti];
                if t.span.is_some() {
                    return Err(bad(row.line, format!("duplicate mRNA ID {id:?}")));
                }
                (t.contig, t.strand) = (row.contig.to_string(), row.strand);
                t.span = Some(row.coords);
                t.attrs = kept_attrs(&row.attrs, &["ID", "Parent"]);
            }
            _ => {
                for parent in attr(&row.attrs, "Parent").unwrap().split(',') {
                    let (gene_id, tx_id) = if let Some(g) = mrna_parent.get(parent) {
                        let gene_id = match g {
                            Some(p) => p.to_string(),
                            None => format!("{parent}.gene"),
                        };
                        (gene_id, parent.to_string())
                    } else if gene_ids.contains(parent) {
                        (parent.to_string(), format!("{parent}.mRNA"))
                    } else if other_ids.contains(parent) {
                        continue; // a child of an unread feature type (tRNA, ncRNA, ...)
                    } else {
                        (format!("{parent}.gene"), parent.to_string())
                    };
                    let gi = gene_of(&mut genes, &mut gene_idx, &gene_id, row);
                    let ti = tx_of(&mut genes, &mut txs, &mut tx_idx, gi, &tx_id, row);
                    let t = &mut txs[ti];
                    match row.kind {
                        Kind::Exon => t.exons.push(row.coords),
                        Kind::Cds => t.cds.push((row.coords, row.phase)),
                        Kind::FiveUtr => t.five.push(row.coords),
                        Kind::ThreeUtr => t.three.push(row.coords),
                        Kind::Gene | Kind::Mrna => unreachable!(),
                    }
                }
            }
        }
    }

    // Finish: sort features, synthesize exons, derive spans.
    let mut txs: Vec<Option<TxAcc>> = txs.into_iter().map(Some).collect();
    let mut out = Vec::with_capacity(genes.len());
    for g in genes {
        let mut transcripts = Vec::with_capacity(g.tx.len());
        for ti in g.tx {
            let mut t = txs[ti].take().unwrap();
            t.cds.sort_by_key(|(c, _)| (c.lend, c.rend));
            t.exons.sort_by_key(|c| (c.lend, c.rend));
            t.five.sort_by_key(|c| (c.lend, c.rend));
            t.three.sort_by_key(|c| (c.lend, c.rend));
            let cds_start_phase = match t.strand {
                Strand::Minus => t.cds.last(),
                _ => t.cds.first(),
            }
            .map_or(0, |&(_, p)| p);
            let cds: Vec<Coordset> = t.cds.iter().map(|&(c, _)| c).collect();
            let exons = if !t.exons.is_empty() {
                t.exons
            } else if !(cds.is_empty() && t.five.is_empty() && t.three.is_empty()) {
                let mut all = cds.clone();
                all.extend_from_slice(&t.five);
                all.extend_from_slice(&t.three);
                merge_coords(all)
            } else {
                vec![t.span.expect("a transcript with no children has a row")]
            };
            transcripts.push(OutTranscript {
                transcript_id: t.id,
                contig: t.contig,
                strand: t.strand,
                exons,
                cds,
                five_utr: t.five,
                three_utr: t.three,
                attrs: t.attrs,
                cds_start_phase,
            });
        }
        let (lend, rend) = match g.span {
            Some(s) => (s.lend, s.rend),
            None => (
                transcripts
                    .iter()
                    .flat_map(|t| &t.exons)
                    .map(|c| c.lend)
                    .min()
                    .unwrap(),
                transcripts
                    .iter()
                    .flat_map(|t| &t.exons)
                    .map(|c| c.rend)
                    .max()
                    .unwrap(),
            ),
        };
        let (contig, strand) = match (g.span, transcripts.first()) {
            (None, Some(t)) => (t.contig.clone(), t.strand),
            _ => (g.contig, g.strand),
        };
        out.push(OutGene {
            gene_id: g.id,
            contig,
            strand,
            lend,
            rend,
            attrs: g.attrs,
            transcripts,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::writer_gff3;

    fn cs(l: i64, r: i64) -> Coordset {
        Coordset::new(l, r)
    }

    fn parse(text: &str) -> Vec<OutGene> {
        parse_gene_models(text, "t.gff3").unwrap()
    }

    fn get<'a>(attrs: &'a Attrs, k: &str) -> Option<&'a Vec<String>> {
        attrs.iter().find(|(key, _)| key == k).map(|(_, v)| v)
    }

    const FULL: &str = "##gff-version 3
chr1\tsrc\tgene\t1\t200\t.\t+\t.\tID=g1;Name=alpha;Note=a%3Bb,c%2Cd
chr1\tsrc\tmRNA\t1\t200\t.\t+\t.\tID=g1.t1;Parent=g1;Dbxref=X:1,Y:2;note=p%3Dq
chr1\tsrc\tfive_prime_UTR\t1\t9\t.\t+\t.\tParent=g1.t1
chr1\tsrc\texon\t101\t200\t.\t+\t.\tParent=g1.t1
chr1\tsrc\texon\t1\t50\t.\t+\t.\tParent=g1.t1
chr1\tsrc\tCDS\t101\t180\t.\t+\t2\tParent=g1.t1
chr1\tsrc\tCDS\t10\t50\t.\t+\t0\tParent=g1.t1
chr1\tsrc\tthree_prime_UTR\t181\t200\t.\t+\t.\tParent=g1.t1
chr1\tsrc\tmRNA\t1\t200\t.\t+\t.\tID=g1.t2;Parent=g1
chr1\tsrc\texon\t1\t200\t.\t+\t.\tParent=g1.t2
chr1\tsrc\tCDS\t10\t180\t.\t+\t0\tParent=g1.t2
";

    #[test]
    fn full_hierarchy_and_attributes() {
        let genes = parse(FULL);
        assert_eq!(genes.len(), 1);
        let g = &genes[0];
        assert_eq!((g.gene_id.as_str(), g.lend, g.rend), ("g1", 1, 200));
        assert_eq!(g.strand, Strand::Plus);
        // gene attrs: file order, ID dropped, multi-values split, escapes decoded
        assert_eq!(g.attrs[0].0, "Name");
        assert_eq!(
            get(&g.attrs, "Note").unwrap(),
            &vec!["a;b".to_string(), "c,d".to_string()]
        );
        assert_eq!(g.transcripts.len(), 2);
        let t = &g.transcripts[0];
        assert_eq!(t.transcript_id, "g1.t1");
        assert_eq!(t.exons, vec![cs(1, 50), cs(101, 200)]); // lend-sorted
        assert_eq!(t.cds, vec![cs(10, 50), cs(101, 180)]);
        assert_eq!(t.five_utr, vec![cs(1, 9)]);
        assert_eq!(t.three_utr, vec![cs(181, 200)]);
        assert_eq!(t.cds_start_phase, 0);
        let keys: Vec<&str> = t.attrs.iter().map(|(k, _)| k.as_str()).collect();
        assert_eq!(keys, vec!["Dbxref", "note"]); // ID/Parent dropped, order kept
        assert_eq!(
            get(&t.attrs, "Dbxref").unwrap(),
            &vec!["X:1".to_string(), "Y:2".to_string()]
        );
        assert_eq!(get(&t.attrs, "note").unwrap(), &vec!["p=q".to_string()]);
        assert_eq!(g.transcripts[1].transcript_id, "g1.t2");
    }

    #[test]
    fn attributes_round_trip_through_the_writer() {
        let genes = parse(FULL);
        let mut buf = Vec::new();
        writer_gff3::write(&mut buf, &genes).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("\tID=g1;Name=alpha;Note=a%3Bb,c%2Cd\n"));
        assert!(out.contains("\tID=g1.t1;Parent=g1;Dbxref=X:1,Y:2;note=p%3Dq\n"));
        // and re-parsing the writer's output is a fixed point
        let again = parse(&out);
        let mut buf2 = Vec::new();
        writer_gff3::write(&mut buf2, &again).unwrap();
        assert_eq!(out, String::from_utf8(buf2).unwrap());
    }

    #[test]
    fn cds_only_mrna_gets_exons_from_cds_and_utrs() {
        let genes = parse(
            "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=m;Parent=g
chr1\ts\tfive_prime_UTR\t1\t9\t.\t+\t.\tParent=m
chr1\ts\tCDS\t10\t40\t.\t+\t0\tParent=m
chr1\ts\tCDS\t60\t90\t.\t+\t1\tParent=m
",
        );
        let t = &genes[0].transcripts[0];
        assert_eq!(t.exons, vec![cs(1, 40), cs(60, 90)]); // UTR+CDS coalesced
        assert_eq!(t.cds, vec![cs(10, 40), cs(60, 90)]);
    }

    #[test]
    fn missing_gene_and_mrna_rows_are_synthesized() {
        let genes = parse(
            "chr1\ts\tmRNA\t10\t90\t.\t-\t.\tID=m1;Parent=gX
chr1\ts\tCDS\t10\t90\t.\t-\t0\tParent=m1
chr1\ts\tmRNA\t200\t300\t.\t+\t.\tID=m2
chr1\ts\tCDS\t200\t300\t.\t+\t0\tParent=m2
chr2\ts\tCDS\t5\t50\t.\t-\t0\tParent=orphan
",
        );
        let ids: Vec<(&str, &str)> = genes
            .iter()
            .map(|g| (g.gene_id.as_str(), g.transcripts[0].transcript_id.as_str()))
            .collect();
        assert_eq!(
            ids,
            vec![("gX", "m1"), ("m2.gene", "m2"), ("orphan.gene", "orphan")]
        );
        // synthesized genes take contig/strand/span from their transcripts
        assert_eq!(
            (genes[0].strand, genes[0].lend, genes[0].rend),
            (Strand::Minus, 10, 90)
        );
        assert_eq!(genes[2].contig, "chr2");
        assert_eq!(genes[2].strand, Strand::Minus);
        assert_eq!(genes[2].transcripts[0].exons, vec![cs(5, 50)]);
    }

    #[test]
    fn cds_parented_by_a_gene_gets_a_synthesized_mrna() {
        let genes = parse(
            "chr1\ts\tgene\t1\t99\t.\t+\t.\tID=g
chr1\ts\tCDS\t1\t99\t.\t+\t0\tParent=g
",
        );
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].transcripts[0].transcript_id, "g.mRNA");
    }

    #[test]
    fn minus_strand_start_phase_comes_from_the_rightmost_cds() {
        let genes = parse(
            "chr1\ts\tgene\t1\t100\t.\t-\t.\tID=g
chr1\ts\tmRNA\t1\t100\t.\t-\t.\tID=m;Parent=g
chr1\ts\tCDS\t60\t100\t.\t-\t2\tParent=m
chr1\ts\tCDS\t1\t40\t.\t-\t1\tParent=m
",
        );
        let t = &genes[0].transcripts[0];
        assert_eq!(t.strand, Strand::Minus);
        assert_eq!(t.cds, vec![cs(1, 40), cs(60, 100)]);
        assert_eq!(t.cds_start_phase, 2);
    }

    #[test]
    fn plus_strand_nonzero_start_phase() {
        let genes = parse(
            "chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=m
chr1\ts\tCDS\t60\t100\t.\t+\t0\tParent=m
chr1\ts\tCDS\t1\t40\t.\t+\t1\tParent=m
",
        );
        assert_eq!(genes[0].transcripts[0].cds_start_phase, 1);
    }

    #[test]
    fn genes_in_first_appearance_order_and_other_types_ignored() {
        let genes = parse(
            "chr2\ts\tgene\t1\t100\t.\t+\t.\tID=b
chr1\ts\tgene\t1\t100\t.\t+\t.\tID=a
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=a1;Parent=a
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=a1
chr2\ts\ttRNA\t1\t100\t.\t+\t.\tID=b1;Parent=b
chr2\ts\texon\t1\t100\t.\t+\t.\tParent=b1
chr2\ts\tmRNA\t1\t100\t.\t+\t.\tID=b2;Parent=b
chr2\ts\texon\t1\t100\t.\t+\t.\tParent=b2
",
        );
        let ids: Vec<&str> = genes.iter().map(|g| g.gene_id.as_str()).collect();
        assert_eq!(ids, vec!["b", "a"]);
        let b: Vec<&str> = genes[0]
            .transcripts
            .iter()
            .map(|t| t.transcript_id.as_str())
            .collect();
        assert_eq!(b, vec!["b2"]); // the tRNA and its exon are not read
    }

    #[test]
    fn bad_rows_are_parse_errors_with_line_numbers() {
        let err = parse_gene_models(
            "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g\nchr1\ts\tmRNA\tX\t100\t.\t+\t.\tID=m\n",
            "f.gff3",
        )
        .map(|_| ())
        .unwrap_err();
        assert!(matches!(err, CombinrError::Parse { line: 2, .. }), "{err}");
        let err = parse_gene_models("chr1\ts\tCDS\t1\t9\t.\t*\t0\tParent=m\n", "f")
            .map(|_| ())
            .unwrap_err();
        assert!(err.to_string().contains("bad strand"), "{err}");
        let err = parse_gene_models("chr1\ts\tCDS\t1\t9\t.\t+\t3\tParent=m\n", "f")
            .map(|_| ())
            .unwrap_err();
        assert!(err.to_string().contains("bad phase"), "{err}");
    }
}
