//! Expand a sequence of CIGAR operations using the reference and the string to split
//! match elements into sequence match and sequence mismatch elements.
//!
//! Expand a CIGAR string, using the reference and the sequence to split
//! match elements into sequence match and sequence mismatch elements.
//! # Example
//!
//! ```rust
//! use cigar_utils::{CigarOp, CigarElement};
//! use cigar_utils::expand::expand_cigar_operations;
//! 
//! let reference = b"ACGT";
//! let reference = b"ACGT";
//! let seq = b"AGGT";
//! let cigar = "4M";
//!
//! let expanded = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
//! let cigar_string = CigarElement::cigar_string(expanded);
//!
//! // The expanded result will split the match into equal and diff elements:
//! assert_eq!(cigar_string, "1=1X2=");
//! ```

use crate::{CigarElement, CigarIterator, CigarOp, error::CigarError};

/// Expand a CIGAR string, using the reference and the sequence to split
/// match elements into sequence match and sequence mismatch elements.
pub fn expand_cigar_operations<R: AsRef<[u8]>, S: AsRef<[u8]>>(
    reference_position: usize,
    cigar: &str,
    reference: &R,
    seq: &S,
) -> std::result::Result<Vec<CigarElement>, CigarError> {
    let mut expanded = Vec::new();
    let mut reference_position = reference_position;
    let mut read_sequence_position = 0;

    for elem in CigarIterator::new(cigar) {
        let elem = elem?;
        match elem.op {
            CigarOp::Match => {
                // Split the match element into sequence match and mismatch elements
                let seq_slice = &seq.as_ref()[read_sequence_position..read_sequence_position + elem.length as usize];
                let ref_slice = &reference.as_ref()[reference_position..reference_position + elem.length as usize];
                let mut match_length = 0;
                let mut mismatch_length = 0;
                for (s, r) in seq_slice.iter().zip(ref_slice.iter()) {
                    if s == r {
                        if mismatch_length > 0 {
                            expanded.push(CigarElement::new(mismatch_length, CigarOp::Diff));
                            mismatch_length = 0;
                        }
                        match_length += 1;
                    } else {
                        if match_length > 0 {
                            expanded.push(CigarElement::new(match_length, CigarOp::Equal));
                            match_length = 0;
                        }
                        mismatch_length += 1;
                    }
                }
                if match_length > 0 {
                    expanded.push(CigarElement::new(match_length, CigarOp::Equal));
                }
                if mismatch_length > 0 {
                    expanded.push(CigarElement::new(mismatch_length, CigarOp::Diff));
                }
                read_sequence_position += elem.length as usize;
                reference_position += elem.length as usize;
            }
            CigarOp::Insertion => {
                read_sequence_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::Deletion => {
                reference_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::Skip => {
                reference_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::SoftClip => {
                read_sequence_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::HardClip => {
                // We want the read sequence position only.
                // read_sequence_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::Padding => {
                reference_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::Equal => {
                read_sequence_position += elem.length as usize;
                reference_position += elem.length as usize;
                expanded.push(elem);
            },
            CigarOp::Diff => {
                read_sequence_position += elem.length as usize;
                reference_position += elem.length as usize;
                expanded.push(elem);
            },
        }
    }

    Ok(expanded)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CigarOp;

    #[test]
    fn test_expand_cigar_all_match() {
        let reference = b"ACGT";
        let seq = b"ACGT";
        let cigar = "4M";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].op, CigarOp::Equal);
        assert_eq!(result[0].length, 4);
    }

    #[test]
    fn test_expand_cigar_all_mismatch() {
        let reference = b"ACGT";
        let seq = b"TGCA";
        let cigar = "4M";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].op, CigarOp::Diff);
        assert_eq!(result[0].length, 4);
    }

    #[test]
    fn test_expand_cigar_mixed_match_mismatch() {
        let reference = b"ACGT";
        let seq = b"AGGT";
        let cigar = "4M";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        // A==A, C!=G, G==G, T==T
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].op, CigarOp::Equal);
        assert_eq!(result[0].length, 1);
        assert_eq!(result[1].op, CigarOp::Diff);
        assert_eq!(result[1].length, 1);
        assert_eq!(result[2].op, CigarOp::Equal);
        assert_eq!(result[2].length, 2);
    }

    #[test]
    fn test_expand_cigar_with_insertion_and_deletion() {
        let reference = b"ACGT";
        let seq = b"ACGTT";
        let cigar = "4M1I";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        // 4M should be all match, 1I is insertion
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].op, CigarOp::Equal);
        assert_eq!(result[0].length, 4);
        assert_eq!(result[1].op, CigarOp::Insertion);
        assert_eq!(result[1].length, 1);
    }

    #[test]
    fn test_expand_cigar_with_softclip() {
        let reference = b"ACGT";
        let seq = b"AACGT";
        let cigar = "1S4M";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].op, CigarOp::SoftClip);
        assert_eq!(result[0].length, 1);
        assert_eq!(result[1].op, CigarOp::Equal);
        assert_eq!(result[1].length, 4);
    }

    #[test]
    fn test_expand_cigar_with_right_softclip() {
        let reference = b"ACGT";
        let seq = b"ACGTT";
        let cigar = "4M1S";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].op, CigarOp::Equal);
        assert_eq!(result[0].length, 4);
        assert_eq!(result[1].op, CigarOp::SoftClip);
        assert_eq!(result[1].length, 1);
    }

    #[test]
    fn test_expand_cigar_with_hardclip() {
        let reference = b"ACGT";
        let seq = b"ACGT"; // Hard clipped base is not present in sequence
        let cigar = "1H4M";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].op, CigarOp::HardClip);
        assert_eq!(result[0].length, 1);
        assert_eq!(result[1].op, CigarOp::Equal);
        assert_eq!(result[1].length, 4);
    }

    #[test]
    fn test_expand_cigar_with_right_hardclip() {
        let reference = b"ACGT";
        let seq = b"ACGT"; // Only first base present, rest hard clipped
        let cigar = "4M1H";
        let result = expand_cigar_operations(0, cigar, &reference, &seq).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].op, CigarOp::Equal);
        assert_eq!(result[0].length, 4);
        assert_eq!(result[1].op, CigarOp::HardClip);
        assert_eq!(result[1].length, 1);
    }
}
