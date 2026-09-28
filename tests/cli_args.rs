//! CLI argument validation for `assemble`'s --gene-pred / --models / --genome options.

use assert_cmd::Command;
use predicates::prelude::*;
use predicates::str::contains;

#[test]
fn lone_genome_errors_and_names_both_anchors() {
    Command::cargo_bin("combinr")
        .unwrap()
        .args(["assemble", "-i", "x.gff3", "--genome", "g.fa"])
        .assert()
        .failure()
        .stderr(contains("--gene-pred").and(contains("--models")));
}

#[test]
fn models_conflicts_with_gene_pred() {
    Command::cargo_bin("combinr")
        .unwrap()
        .args([
            "assemble",
            "-i",
            "x.gff3",
            "--genome",
            "g.fa",
            "--models",
            "m.gff3",
            "--gene-pred",
            "p.gff3",
        ])
        .assert()
        .failure()
        .stderr(contains("cannot be used with"));
}

#[test]
fn models_requires_genome() {
    Command::cargo_bin("combinr")
        .unwrap()
        .args(["assemble", "-i", "x.gff3", "--models", "m.gff3"])
        .assert()
        .failure()
        .stderr(contains("--genome"));
}
