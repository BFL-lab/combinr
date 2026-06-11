//! Golden parity for the assembler core (M1): for every
//! `pasa_cpp_sample_input*`, run both the compiled C++ `pasa` binary and
//! `combinr assemble-tokens` and assert identical assembly *sets* (member sets +
//! merged structures). Sets, not line order, because tied `lend` values make the
//! C++ `std::sort` (unstable) and Rust's stable sort number members differently
//! without ever changing which assemblies are produced.
//!
//! Skips with a notice if the reference `pasa` binary cannot be found (build it
//! with `make` in `pasa_cpp`, or point `$COMBINR_PASA_DIR` at the PASA root).

mod common;

use std::process::Command;

#[test]
fn assembler_matches_golden_pasa_on_all_samples() {
    let Some(pasa) = common::pasa_binary() else {
        eprintln!("SKIP: reference `pasa` binary not found (set $COMBINR_PASA_DIR or build it)");
        return;
    };
    let dir = pasa.parent().unwrap().to_path_buf();
    let combinr = env!("CARGO_BIN_EXE_combinr");

    let inputs = common::sample_inputs(&dir);
    assert!(!inputs.is_empty(), "no sample inputs in {}", dir.display());

    let mut checked = 0;
    for input in &inputs {
        let golden = Command::new(&pasa).arg(input).output().expect("run pasa");
        let expected = common::parse_assemblies(&String::from_utf8_lossy(&golden.stdout));

        let mine = Command::new(combinr)
            .args(["assemble-tokens", input.to_str().unwrap()])
            .output()
            .expect("run combinr");
        assert!(
            mine.status.success(),
            "combinr failed on {}: {}",
            input.display(),
            String::from_utf8_lossy(&mine.stderr)
        );
        let got = common::parse_assemblies(&String::from_utf8_lossy(&mine.stdout));

        assert_eq!(
            got,
            expected,
            "assembly set mismatch on {}\n  only in combinr: {:?}\n  only in pasa: {:?}",
            input.display(),
            got.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&got).collect::<Vec<_>>(),
        );
        checked += 1;
    }
    assert!(checked > 0);
    eprintln!("golden parity: {checked} sample inputs matched");
}
