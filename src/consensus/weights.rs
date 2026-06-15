//! EVM weights file + evidence classes.
//!
//! The weights file is three whitespace-separated columns: `EV_CLASS  EV_TYPE  WEIGHT`,
//! where `EV_TYPE` matches GFF column 2 (the "source") of the corresponding evidence
//! rows. Port of `evidence_modeler.pl::readWeights` (~line 2883). The four classes and
//! their scoring roles match `%ALLOWABLE_EVIDENCE_CLASSES` (~line 87).

use crate::error::{CombinrError, Result};
use std::collections::HashMap;
use std::path::Path;

/// EVM evidence class. Governs how an evidence type contributes to scoring.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub enum EvClass {
    Protein,
    Transcript,
    AbinitioPrediction,
    OtherPrediction,
}

impl EvClass {
    pub fn parse(s: &str) -> Option<EvClass> {
        match s {
            "PROTEIN" => Some(EvClass::Protein),
            "TRANSCRIPT" => Some(EvClass::Transcript),
            "ABINITIO_PREDICTION" => Some(EvClass::AbinitioPrediction),
            "OTHER_PREDICTION" => Some(EvClass::OtherPrediction),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            EvClass::Protein => "PROTEIN",
            EvClass::Transcript => "TRANSCRIPT",
            EvClass::AbinitioPrediction => "ABINITIO_PREDICTION",
            EvClass::OtherPrediction => "OTHER_PREDICTION",
        }
    }

    /// PROTEIN + ABINITIO paint the per-base coding vector
    /// (`add_match_coverage` gate, evidence_modeler.pl ~line 2602).
    pub fn paints_coding_vector(self) -> bool {
        matches!(self, EvClass::Protein | EvClass::AbinitioPrediction)
    }

    /// Only ABINITIO contributes to intergenic scoring
    /// (`%PREDICTION_PROGS_CONTRIBUTE_INTERGENIC`, ~line 199).
    pub fn contributes_intergenic(self) -> bool {
        matches!(self, EvClass::AbinitioPrediction)
    }

    /// TRANSCRIPT + OTHER_PREDICTION contribute exon-specific exact-match scores
    /// instead of coding-vector coverage (`score_exons`, ~line 2238).
    pub fn exon_specific(self) -> bool {
        matches!(self, EvClass::Transcript | EvClass::OtherPrediction)
    }

    /// ABINITIO + OTHER_PREDICTION provide trusted initial/terminal exon types;
    /// PROTEIN/TRANSCRIPT only seed internal exons + boundaries (~line 90).
    pub fn provides_terminal_exon_types(self) -> bool {
        matches!(self, EvClass::AbinitioPrediction | EvClass::OtherPrediction)
    }

    /// A gene-prediction kind (parsed from CDS rows, with an authoritative frame).
    pub fn is_prediction(self) -> bool {
        matches!(self, EvClass::AbinitioPrediction | EvClass::OtherPrediction)
    }
}

/// Evidence weights, keyed by `ev_type` (== GFF column-2 source string).
#[derive(Clone, Debug, Default)]
pub struct Weights {
    by_type: HashMap<String, (EvClass, f64)>,
}

impl Weights {
    pub fn parse_file(path: &Path) -> Result<Weights> {
        let text = std::fs::read_to_string(path).map_err(|e| CombinrError::Parse {
            file: path.display().to_string(),
            line: 0,
            msg: format!("cannot read weights file: {e}"),
        })?;
        Self::parse_str(&text, &path.display().to_string())
    }

    pub fn parse_str(text: &str, file: &str) -> Result<Weights> {
        let mut by_type = HashMap::new();
        for (lineno, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 3 {
                return Err(CombinrError::Parse {
                    file: file.to_string(),
                    line: lineno + 1,
                    msg: format!(
                        "expected 3 columns (CLASS TYPE WEIGHT), got {}",
                        fields.len()
                    ),
                });
            }
            let class = EvClass::parse(fields[0]).ok_or_else(|| CombinrError::Parse {
                file: file.to_string(),
                line: lineno + 1,
                msg: format!("unknown evidence class {:?}", fields[0]),
            })?;
            let weight: f64 = fields[2].parse().map_err(|_| CombinrError::Parse {
                file: file.to_string(),
                line: lineno + 1,
                msg: format!("bad weight {:?}", fields[2]),
            })?;
            by_type.insert(fields[1].to_string(), (class, weight));
        }
        Ok(Weights { by_type })
    }

    /// `(class, weight)` for an evidence type, if present.
    pub fn lookup(&self, ev_type: &str) -> Option<(EvClass, f64)> {
        self.by_type.get(ev_type).copied()
    }

    pub fn class(&self, ev_type: &str) -> Option<EvClass> {
        self.by_type.get(ev_type).map(|&(c, _)| c)
    }

    pub fn weight(&self, ev_type: &str) -> Option<f64> {
        self.by_type.get(ev_type).map(|&(_, w)| w)
    }

    pub fn len(&self) -> usize {
        self.by_type.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_type.is_empty()
    }

    /// Sum of all gene-prediction weights (EVM's `SUM_GENEPRED_WEIGHTS`, used later
    /// for start/stop peak-augmentation thresholds).
    pub fn sum_genepred_weights(&self) -> f64 {
        self.by_type
            .values()
            .filter(|&&(c, _)| c.is_prediction())
            .map(|&(_, w)| w)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
# example weights
PROTEIN nap 1
PROTEIN genewise 5
TRANSCRIPT gap2 1
TRANSCRIPT alignAssembly 10
ABINITIO_PREDICTION fgenesh 1
ABINITIO_PREDICTION genemark 1
OTHER_PREDICTION ensembl 6
";

    #[test]
    fn parses_sample() {
        let w = Weights::parse_str(SAMPLE, "weights.txt").unwrap();
        assert_eq!(w.len(), 7);
        assert_eq!(w.lookup("genewise"), Some((EvClass::Protein, 5.0)));
        assert_eq!(w.class("ensembl"), Some(EvClass::OtherPrediction));
        assert_eq!(w.weight("alignAssembly"), Some(10.0));
        assert_eq!(w.lookup("missing"), None);
    }

    #[test]
    fn sum_genepred_weights_counts_predictions_only() {
        // abinitio fgenesh(1) + genemark(1) + other ensembl(6) = 8
        let w = Weights::parse_str(SAMPLE, "w").unwrap();
        assert_eq!(w.sum_genepred_weights(), 8.0);
    }

    #[test]
    fn unknown_class_errors() {
        assert!(Weights::parse_str("BOGUS x 1", "w").is_err());
    }

    #[test]
    fn too_few_columns_errors() {
        assert!(Weights::parse_str("PROTEIN nap", "w").is_err());
    }

    #[test]
    fn class_helpers() {
        assert!(EvClass::Protein.paints_coding_vector());
        assert!(EvClass::AbinitioPrediction.paints_coding_vector());
        assert!(!EvClass::Transcript.paints_coding_vector());
        assert!(EvClass::AbinitioPrediction.contributes_intergenic());
        assert!(!EvClass::OtherPrediction.contributes_intergenic());
        assert!(EvClass::Transcript.exon_specific());
        assert!(EvClass::OtherPrediction.provides_terminal_exon_types());
        assert_eq!(
            EvClass::parse(EvClass::Protein.as_str()),
            Some(EvClass::Protein)
        );
    }
}
