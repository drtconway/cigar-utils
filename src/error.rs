//! CIGAR errors.

use std::error::Error;
use std::fmt::Display;

/// Errors that can occur while parsing CIGAR strings.
#[derive(Debug)]
pub enum CigarError {
    /// An error indicating an invalid character in the CIGAR string.
    InvalidCharacter(char),
    /// An error indicating a missing count in a CIGAR element.
    MissingCount(char),
    /// An error indicating a missing operation in a CIGAR element.
    MissingOperation(u32),
    /// An external error.
    External(Box<dyn Error + Send + Sync + 'static>),
}

impl Display for CigarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CigarError::InvalidCharacter(c) => write!(f, "Invalid character in CIGAR string: {}", c),
            CigarError::MissingCount(c) => write!(f, "Missing count in CIGAR element (found '{}')", c),
            CigarError::MissingOperation(length) => write!(f, "Missing operation in CIGAR element (length was {})", length),
            CigarError::External(_) => write!(f, "External error"),
        }
    }
}

impl Error for CigarError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            CigarError::External(e) => Some(&**e),
            _ => None,
        }
    }
}
