/// Global program name
///
pub const PROGRAM_NAME: &str = env!("CARGO_PKG_NAME");

/// Global version number
///
/// All client code should refer directly to this copy instead of using various possibly conflicting environment variables
pub const PROGRAM_VERSION: &str = env!("VERGEN_GIT_DESCRIBE");
