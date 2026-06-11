//! Alt-splice golden parity (Algorithm 2) against committed PASA references.
//!
//! For each vendored input under `tests/data/altsplice/`, PASA's real
//! `CDNA::Alternative_splice_comparer` was run on combinr's emitted isoforms and
//! its events committed as `*.events.golden`. combinr's own classification must
//! reproduce that exact event set. No PASA code runs here (the regeneration
//! harness lives in `validation/`).

mod common;

use combinr::filter::Filters;
use combinr::pipeline::analyze_sources;

#[test]
fn altsplice_matches_pasa_golden() {
    let fixtures = common::altsplice_fixtures();
    assert!(!fixtures.is_empty(), "no alt-splice fixtures vendored");
    for (input, golden) in fixtures {
        let result = analyze_sources(std::slice::from_ref(&input), 20, &Filters::none()).unwrap();
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
