//! Augmented CIGAR operations.
//! 
//! Augmented CIGAR operations provide additional context to the standard CIGAR operations by including read and reference positions.
//! 
//! This module also provides iterators over sequences of them derived from an alignment position and a cigar string.

use crate::error::CigarError;
use crate::{CigarElement, CigarIterator, CigarOp};

/// An augmented CIGAR operation element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AugmentedCigarElement {
    /// The length of the CIGAR operation.
    pub length: u32,
    /// The type of the CIGAR operation.
    pub op: CigarOp,
    /// The read position of the CIGAR operation.
    pub read_position: u32,
    /// The reference position of the CIGAR operation.
    pub reference_position: u32,
}

impl Ord for AugmentedCigarElement {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.reference_position.cmp(&other.reference_position) {
            std::cmp::Ordering::Equal => {
                match self.op.cmp(&other.op) {
                    std::cmp::Ordering::Equal => self.length.cmp(&other.length),
                    ord => ord,
                }
            }
            ord => ord,
        }
    }
}

impl PartialOrd for AugmentedCigarElement {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// An iterator over augmented CIGAR elements.
pub struct AugmentedCigarIterator<'a> {
    pub(crate) inner: CigarIterator<'a>,
    pub(crate) read_position: u32,
    pub(crate) reference_position: u32,
}

impl<'a> From<(CigarIterator<'a>, u32)> for AugmentedCigarIterator<'a> {
    fn from(value: (CigarIterator<'a>, u32)) -> Self {
        let (inner, reference_position) = value;
        AugmentedCigarIterator {
            inner,
            read_position: 0,
            reference_position,
        }
    }
}

impl<'a> From<(&'a str, u32)> for AugmentedCigarIterator<'a> {
    fn from(value: (&'a str, u32)) -> Self {
        let (cigar_str, reference_position) = value;
        let inner = CigarIterator {
            chars: cigar_str.chars(),
        };
        AugmentedCigarIterator {
            inner,
            read_position: 0,
            reference_position,
        }
    }
}

impl<'a> Iterator for AugmentedCigarIterator<'a> {
    type Item = std::result::Result<AugmentedCigarElement, CigarError>;

    fn next(&mut self) -> Option<Self::Item> {
        let inner_elem = self.inner.next()?;
        match inner_elem {
            Ok(CigarElement { length, op }) => {
                let read_position = self.read_position;
                let reference_position = self.reference_position;
                let elem = AugmentedCigarElement {
                    length,
                    op,
                    read_position,
                    reference_position,
                };
                match op {
                    CigarOp::Match => {
                        self.read_position += length;
                        self.reference_position += length;
                    }
                    CigarOp::Insertion => {
                        self.read_position += length;
                    }
                    CigarOp::Deletion => {
                        self.reference_position += length;
                    }
                    CigarOp::Skip => {
                        self.reference_position += length;
                    }
                    CigarOp::SoftClip => {
                        self.read_position += length;
                    }
                    CigarOp::HardClip => {
                        self.read_position += length;
                    }
                    CigarOp::Padding => {
                        self.read_position += length;
                    }
                    CigarOp::Equal => {
                        self.read_position += length;
                        self.reference_position += length;
                    }
                    CigarOp::Diff => {
                        self.read_position += length;
                        self.reference_position += length;
                    }
                }
                Some(Ok(elem))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_augmented_cigar_iterator_basic() {
        let cigar = "3M2I4D";
        let iter = AugmentedCigarIterator::from((
            CigarIterator {
                chars: cigar.chars(),
            },
            100,
        ));
        let elems: Vec<_> = iter.collect();
        assert_eq!(elems.len(), 3);
        // 3M
        assert!(matches!(elems[0], Ok(ref e)
        if e.length == 3 && e.op == CigarOp::Match && e.read_position == 0 && e.reference_position == 100));
        // 2I
        assert!(matches!(elems[1], Ok(ref e)
        if e.length == 2 && e.op == CigarOp::Insertion && e.read_position == 3 && e.reference_position == 103));
        // 4D
        assert!(matches!(elems[2], Ok(ref e)
        if e.length == 4 && e.op == CigarOp::Deletion && e.read_position == 5 && e.reference_position == 103));
    }

    #[test]
    fn test_augmented_cigar_iterator_positions() {
        let cigar = "2M1I2D1M";
        let iter = AugmentedCigarIterator::from((
            CigarIterator {
                chars: cigar.chars(),
            },
            50,
        ));
        let elems: Vec<_> = iter.collect();
        assert_eq!(elems.len(), 4);
        assert!(matches!(elems[0], Ok(ref e)
        if e.length == 2 && e.op == CigarOp::Match && e.read_position == 0 && e.reference_position == 50));
        assert!(matches!(elems[1], Ok(ref e)
        if e.length == 1 && e.op == CigarOp::Insertion && e.read_position == 2 && e.reference_position == 52));
        assert!(matches!(elems[2], Ok(ref e)
        if e.length == 2 && e.op == CigarOp::Deletion && e.read_position == 3 && e.reference_position == 52));
        assert!(matches!(elems[3], Ok(ref e)
        if e.length == 1 && e.op == CigarOp::Match && e.read_position == 3 && e.reference_position == 54));
    }

    #[test]
    fn test_augmented_cigar_iterator_error_propagation() {
        let cigar = "2M1Z";
        let iter = AugmentedCigarIterator::from((
            CigarIterator {
                chars: cigar.chars(),
            },
            0,
        ));
        let elems: Vec<_> = iter.collect();
        assert_eq!(elems.len(), 2);
        assert!(matches!(elems[0], Ok(ref e) if e.length == 2 && e.op == CigarOp::Match));
        assert!(matches!(elems[1], Err(CigarError::InvalidCharacter('Z'))));
    }

    #[test]
    fn test_augmented_cigar_iterator_from_str() {
        let cigar = "1M2I";
        let iter = AugmentedCigarIterator::from((cigar, 10));
        let elems: Vec<_> = iter.collect();
        assert_eq!(elems.len(), 2);
        assert!(matches!(elems[0], Ok(ref e)
        if e.length == 1 && e.op == CigarOp::Match && e.read_position == 0 && e.reference_position == 10));
        assert!(matches!(elems[1], Ok(ref e)
        if e.length == 2 && e.op == CigarOp::Insertion && e.read_position == 1 && e.reference_position == 11));
    }
}
