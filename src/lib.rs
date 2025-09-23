//! A Rust library for parsing and working with CIGAR strings in bioinformatics.
//! 
//! This crate provides types and utilities for handling CIGAR operations, including parsing, iteration, and error handling.
//! 
//! # Features
//! - Iterator for parsing CIGAR strings
//! - Augmented CIGAR operations that contextualize the individual operations to an alignment.
//! - Collation of multiple augmented CIGAR operations across multiple CIGAR strings.

#![deny(missing_docs)]

use std::convert::TryFrom;
use std::fmt::Display;

pub mod augmented_cigar;
pub mod collated;
pub mod error;
pub mod expand;

/// CIGAR operation types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CigarOp {
    /// Alignment match (can be a sequence match or mismatch) (M).
    Match,
    /// Insertion to the reference (I).
    Insertion,
    /// Deletion from the reference (D).
    Deletion,
    /// Skipped region from the reference (N).
    Skip,
    /// Soft clipping (clipped sequences present in SEQ) (S).
    SoftClip,
    /// Hard clipping (clipped sequences NOT present in SEQ) (H).
    HardClip,
    /// Padding (silent deletion from padded reference) (P).
    Padding,
    /// Sequence match (=).
    Equal,
    /// Sequence mismatch (X).
    Diff,
}

impl Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let c = match self {
            CigarOp::Match => 'M',
            CigarOp::Insertion => 'I',
            CigarOp::Deletion => 'D',
            CigarOp::Skip => 'N',
            CigarOp::SoftClip => 'S',
            CigarOp::HardClip => 'H',
            CigarOp::Padding => 'P',
            CigarOp::Equal => '=',
            CigarOp::Diff => 'X',
        };
        write!(f, "{}", c)
    }
}

impl From<CigarOp> for u8 {
    fn from(op: CigarOp) -> u8 {
        match op {
            CigarOp::Match => 0,
            CigarOp::Insertion => 1,
            CigarOp::Deletion => 2,
            CigarOp::Skip => 3,
            CigarOp::SoftClip => 4,
            CigarOp::HardClip => 5,
            CigarOp::Padding => 6,
            CigarOp::Equal => 7,
            CigarOp::Diff => 8,
        }
    }
}

impl TryFrom<u8> for CigarOp {
    type Error = u8;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(CigarOp::Match),
            1 => Ok(CigarOp::Insertion),
            2 => Ok(CigarOp::Deletion),
            3 => Ok(CigarOp::Skip),
            4 => Ok(CigarOp::SoftClip),
            5 => Ok(CigarOp::HardClip),
            6 => Ok(CigarOp::Padding),
            7 => Ok(CigarOp::Equal),
            8 => Ok(CigarOp::Diff),
            _ => Err(value),
        }
    }
}

/// A single CIGAR operation element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CigarElement {
    /// The length of the CIGAR operation.
    pub length: u32,
    /// The type of the CIGAR operation.
    pub op: CigarOp,
}

impl CigarElement {
    /// Create a new CIGAR element.
    pub fn new(length: u32, op: CigarOp) -> Self {
        CigarElement { length, op }
    }

    /// Convert a sequence of CIGAR elements into a CIGAR string.
    pub fn cigar_string<V: IntoIterator<Item = CigarElement>>(elements: V) -> String {
        elements.into_iter().map(|e| format!("{}", e)).collect()
    }
}

impl Display for CigarElement {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}{}", self.length, self.op)
    }
}

/// An iterator over CIGAR elements.
pub struct CigarIterator<'a> {
    chars: std::str::Chars<'a>,
}

impl<'a> CigarIterator<'a> {
    /// Create a new CIGAR iterator.
    pub fn new(cigar: &'a str) -> Self {
        CigarIterator { chars: cigar.chars() }
    }
}

impl<'a> Iterator for CigarIterator<'a> {
    type Item = std::result::Result<CigarElement, error::CigarError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut digit_count = 0;
        let mut length = 0;

        while let Some(c) = self.chars.next() {
            match c {
                '0'..='9' => {
                    length = length * 10 + (c as u32 - '0' as u32);
                    digit_count += 1;
                    continue;
                }
                _ => {}
            }
            if digit_count == 0 {
                return Some(Err(error::CigarError::MissingCount(c)));
            }

            match c {
                'M' => return Some(Ok(CigarElement::new(length, CigarOp::Match))),
                'I' => return Some(Ok(CigarElement::new(length, CigarOp::Insertion))),
                'D' => return Some(Ok(CigarElement::new(length, CigarOp::Deletion))),
                'N' => return Some(Ok(CigarElement::new(length, CigarOp::Skip))),
                'S' => return Some(Ok(CigarElement::new(length, CigarOp::SoftClip))),
                'H' => return Some(Ok(CigarElement::new(length, CigarOp::HardClip))),
                'P' => return Some(Ok(CigarElement::new(length, CigarOp::Padding))),
                '=' => return Some(Ok(CigarElement::new(length, CigarOp::Equal))),
                'X' => return Some(Ok(CigarElement::new(length, CigarOp::Diff))),
                _ => {
                    return Some(Err(error::CigarError::InvalidCharacter(c)));
                }
            }
        }

        if digit_count > 0 {
            return Some(Err(error::CigarError::MissingOperation(length)));
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use crate::error::CigarError;

    use super::*;

    #[test]
    fn test_cigar_iterator_basic() {
        let cigar = "10M5I3D";
        let iter = CigarIterator {
            chars: cigar.chars(),
        };
        let elems: Vec<_> = iter.collect();
        assert_eq!(elems.len(), 3);
        assert!(matches!(elems[0], Ok(ref e) if e.length == 10 && matches!(e.op, CigarOp::Match)));
        assert!(
            matches!(elems[1], Ok(ref e) if e.length == 5 && matches!(e.op, CigarOp::Insertion))
        );
        assert!(
            matches!(elems[2], Ok(ref e) if e.length == 3 && matches!(e.op, CigarOp::Deletion))
        );
    }

    #[test]
    fn test_cigar_iterator_all_ops() {
        let cigar = "1M2I3D4N5S6H7P8=9X";
        let iter = CigarIterator {
            chars: cigar.chars(),
        };
        let expected = [
            (1, CigarOp::Match),
            (2, CigarOp::Insertion),
            (3, CigarOp::Deletion),
            (4, CigarOp::Skip),
            (5, CigarOp::SoftClip),
            (6, CigarOp::HardClip),
            (7, CigarOp::Padding),
            (8, CigarOp::Equal),
            (9, CigarOp::Diff),
        ];
        for (i, elem) in iter.enumerate() {
            match elem {
                Ok(e) => {
                    assert_eq!(e.length, expected[i].0);
                    assert!(
                        std::mem::discriminant(&e.op) == std::mem::discriminant(&expected[i].1)
                    );
                }
                Err(_) => panic!("Unexpected error at index {}", i),
            }
        }
    }

    #[test]
    fn test_cigar_iterator_invalid_char() {
        let cigar = "5M2Z";
        let iter = CigarIterator {
            chars: cigar.chars(),
        };
        let elems: Vec<_> = iter.collect();
        assert!(matches!(elems[0], Ok(_)));
        assert!(matches!(elems[1], Err(CigarError::InvalidCharacter('Z'))));
    }

    #[test]
    fn test_cigar_iterator_missing_count() {
        let cigar = "M5I";
        let iter = CigarIterator {
            chars: cigar.chars(),
        };
        let elems: Vec<_> = iter.collect();
        assert!(matches!(elems[0], Err(CigarError::MissingCount('M'))));
        assert!(
            matches!(elems[1], Ok(ref e) if e.length == 5 && matches!(e.op, CigarOp::Insertion))
        );
    }
}
