//! Alt-splice golden parity (Algorithm 2) against committed PASA references.
//!
//! For each vendored input under `tests/data/altsplice/`, PASA's real
//! `CDNA::Alternative_splice_comparer` was run on combinr's emitted isoforms and
//! its events committed as `*.events.golden`. combinr's own classification must
//! reproduce that exact event set. No PASA code runs here; the one-off
//! generation harness is preserved in git history.
//!
//! Scope: because the golden is PASA's comparer applied to *combinr's* isoforms,
//! this validates the alt-splice **classification logic** against PASA — not
//! combinr's isoform/locus construction. The underlying assembly is validated
//! separately against PASA's C++ assembler in `golden_assembler.rs`; locus
//! grouping is combinr-specific and covered by unit tests in `altsplice::locus`.

mod common;

use combinr::filter::Filters;
use combinr::pipeline::analyze_sources;

#[test]
fn altsplice_matches_pasa_golden() {
    let fixtures = common::altsplice_fixtures();
    assert!(!fixtures.is_empty(), "no alt-splice fixtures vendored");
    for (input, golden) in fixtures {
        let result =
            analyze_sources(std::slice::from_ref(&input), 20, 0.0, &Filters::none()).unwrap();
        let got = common::canon_events(&result.events);
        let expected = common::load_lines(&golden);
        assert_eq!(
            got,
            expected,
            "alt-splice event mismatch on {}\n  only in combinr: {:?}\n  only in PASA golden: {:?}",
            input.display(),
            got.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&got).collect::<Vec<_>>(),
        );
    }
}
