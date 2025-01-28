use vergen_gitcl::{Emitter, GitclBuilder};

fn main() {
    let gitcl = GitclBuilder::default()
        .describe(true, true, None)
        .build()
        .unwrap();
    Emitter::default()
        .fail_on_error()
        .add_instructions(&gitcl)
        .unwrap()
        .emit()
        .unwrap();

    // emit build handles the git configuration and build.rs, but we also need to track the toml and src folder to catch dirty
    println!("cargo:rerun-if-changed=Cargo.toml");
    println!("cargo:rerun-if-changed=src");
}
