//! Golden parity for the orientation wrapper (M2), via the library API. The
//! sample inputs are single-orientation clusters, so the wrapper's two forced
//! passes are identical and its greedy cover must reproduce the reference
//! assembly set exactly. Mixed-orientation/flex behavior is covered by unit
//! tests in `assemble::wrapper`.

mod common;

use combinr::assemble::assemble_cluster;
use combinr::model::Strand;
use combinr::token::parse_tokens;
use std::collections::BTreeSet;
use std::process::Command;

#[test]
fn wrapper_matches_golden_pasa_on_all_samples() {
    let Some(pasa) = common::pasa_binary() else {
        eprintln!("SKIP: reference `pasa` binary not found");
        return;
    };
    let dir = pasa.parent().unwrap().to_path_buf();

    for input in common::sample_inputs(&dir) {
        let text = std::fs::read_to_string(&input).unwrap();
        let aligns = parse_tokens(&text, input.to_str().unwrap()).unwrap();
        let asms = assemble_cluster(&aligns, 20).unwrap();

        let got: BTreeSet<common::CanonAssembly> = asms
            .iter()
            .map(|a| {
                let members: BTreeSet<String> = a.contained_accs.iter().cloned().collect();
                let mut fields = vec![match a.orient {
                    Strand::Plus => "+".to_string(),
                    Strand::Minus => "-".to_string(),
                    Strand::Unknown => "?".to_string(),
                }];
                for s in &a.structure.segments {
                    fields.push(format!("{}-{}", s.coords.lend, s.coords.rend));
                }
                (members, fields)
            })
            .collect();

        let golden = Command::new(&pasa).arg(&input).output().expect("run pasa");
        let expected = common::parse_assemblies(&String::from_utf8_lossy(&golden.stdout));

        assert_eq!(
            got,
            expected,
            "wrapper assembly set mismatch on {}\n  only in combinr: {:?}\n  only in pasa: {:?}",
            input.display(),
            got.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&got).collect::<Vec<_>>(),
        );
    }
}
