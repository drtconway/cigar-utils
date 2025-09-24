//! Collated CIGAR operations.
//!
//! This module provides a collated iterator over augmented CIGAR elements, allowing for
//! efficient processing and analysis of CIGAR strings across multiple alignments.
//!
//! # Example: Collating CIGAR Strings
//!
//! ```rust
//! use cigar_utils::collated::CollatedAugmentedCigarIterator;
//!
//! // Example input: a vector of CIGAR strings and their starting reference positions
//! let cigars = vec![
//!     std::io::Result::Ok(("2M1I".to_string(), 1, 100)),
//!     std::io::Result::Ok(("1D2M".to_string(), 1, 102)),
//! ];
//!
//! // Create the collated iterator
//! let mut collated = CollatedAugmentedCigarIterator::new(cigars.into_iter());
//!
//! // Iterate and print collated augmented CIGAR elements
//! while let Some(Ok((elem, count))) = collated.next() {
//!     println!(
//!         "ref_pos: {}, op: {:?}, len: {}, count: {}",
//!         elem.reference_position, elem.op, elem.length, count
//!     );
//! }
//! ```
//!
//! This will print each collated event in order of reference position, with the count of how many times each event occurs at that position.

use std::{cmp::Reverse, collections::BinaryHeap, iter::Peekable};

use crate::augmented_cigar::{AugmentedCigarElement, AugmentedCigarIterator};
use crate::error::CigarError;

/// A collated iterator over augmented CIGAR elements.
pub struct CollatedAugmentedCigarIterator<
    Source: Iterator<Item = std::result::Result<(String, u32, u32), E>>,
    E: std::error::Error + Send + Sync + 'static,
> {
    source: Peekable<Source>,
    queue: BinaryHeap<Reverse<AugmentedCigarElement>>,
}

impl<
    Source: Iterator<Item = std::result::Result<(String, u32, u32), E>>,
    E: std::error::Error + Send + Sync + 'static,
> CollatedAugmentedCigarIterator<Source, E>
{
    /// Create a new collated augmented CIGAR iterator.
    pub fn new(source: Source) -> Self {
        let source = source.peekable();
        let queue = BinaryHeap::new();
        CollatedAugmentedCigarIterator { source, queue }
    }
}

impl<
    Source: Iterator<Item = std::result::Result<(String, u32, u32), E>>,
    E: std::error::Error + Send + Sync + 'static,
> Iterator for CollatedAugmentedCigarIterator<Source, E>
{
    type Item = std::result::Result<(AugmentedCigarElement, usize), CigarError>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(item) = self.source.peek() {
            let item = match item {
                Ok(ord) => ord,
                Err(_) => {
                    let e = self.source.next().unwrap().unwrap_err();
                    return Some(Err(CigarError::External(Box::new(e))));
                }
            };
            let (cigar_str, chrom_id, reference_position) = item;
            let mut augmented_iter =
                AugmentedCigarIterator::from((cigar_str as &str, *chrom_id, *reference_position))
                    .peekable();
            if let Some(Ok(elem)) = augmented_iter.peek() {
                if let Some(Reverse(existing)) = self.queue.peek() {
                    if elem.chrom_id > existing.chrom_id
                        || (elem.chrom_id == existing.chrom_id
                            && elem.reference_position > existing.reference_position)
                    {
                        break;
                    }
                }
            }
            for elem in augmented_iter {
                match elem {
                    Ok(e) => self.queue.push(Reverse(e)),
                    Err(e) => return Some(Err(e)),
                }
            }
            self.source.next();
        }
        if let Some(Reverse(elem)) = self.queue.pop() {
            let mut count = 1;
            while let Some(Reverse(next)) = self.queue.peek() {
                if *next == elem {
                    self.queue.pop();
                    count += 1;
                } else {
                    break;
                }
            }
            Some(Ok((elem, count)))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::CigarOp;

    use super::*;

    #[test]
    fn test_collated_augmented_cigar_iterator_basic() {
        let cigars = vec![
            std::io::Result::Ok(("2M1I".to_string(), 1, 100)),
            std::io::Result::Ok(("1D2M".to_string(), 1, 102)),
        ];
        let mut collated = CollatedAugmentedCigarIterator::new(cigars.into_iter());
        let mut results = Vec::new();
        while let Some(Ok((elem, count))) = collated.next() {
            results.push((elem, count));
        }
        println!("{:?}", results);
        // Should be sorted by reference_position
        let positions: Vec<_> = results.iter().map(|(e, _)| e.reference_position).collect();
        assert!(positions.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(results.len(), 4);
        // Check the actual elements
        assert_eq!(results[0].0.reference_position, 100); // 2M from first
        assert_eq!(results[1].0.reference_position, 102); // 1I from first and 1D from second (same pos)
        assert_eq!(results[1].0.op, CigarOp::Insertion);
        assert_eq!(results[2].0.reference_position, 102);
        assert_eq!(results[2].0.op, CigarOp::Deletion);
        assert_eq!(results[3].0.reference_position, 103); // 2M from second
    }

    #[test]
    fn test_collated_augmented_cigar_iterator_error() {
        let cigars = vec![
            std::io::Result::Ok(("2M1Z".to_string(), 1, 100)), // Invalid op 'Z'
            std::io::Result::Ok(("1M".to_string(), 1, 101)),
        ];
        let mut collated = CollatedAugmentedCigarIterator::new(cigars.into_iter());
        let mut saw_error = false;
        while let Some(res) = collated.next() {
            match res {
                Ok(_) => {}
                Err(CigarError::InvalidCharacter('Z')) => {
                    saw_error = true;
                    break;
                }
                Err(_) => {}
            }
        }
        assert!(saw_error);
    }
    #[test]
    fn test_collated_augmented_cigar_iterator_chrom_id_collation() {
        let cigars = vec![
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
            std::io::Result::Ok(("1M".to_string(), 1, 101)),
            std::io::Result::Ok(("1M".to_string(), 2, 100)),
            std::io::Result::Ok(("1M".to_string(), 2, 101)),
        ];
        let mut collated = CollatedAugmentedCigarIterator::new(cigars.into_iter());
        let mut results = Vec::new();
        while let Some(Ok((elem, count))) = collated.next() {
            results.push((elem.chrom_id, elem.reference_position, count));
        }
        // Should have two chrom_id values at each position
        assert_eq!(results.len(), 4);
        assert_eq!(results[0], (1, 100, 2));
        assert_eq!(results[1], (1, 101, 1));
        assert_eq!(results[2], (2, 100, 1));
        assert_eq!(results[3], (2, 101, 1));
    }

    #[test]
    fn test_collated_augmented_cigar_iterator_chrom_id_grouping() {
        let cigars = vec![
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
            std::io::Result::Ok(("1M".to_string(), 2, 100)),
            std::io::Result::Ok(("1M".to_string(), 2, 100)),
        ];
        let mut collated = CollatedAugmentedCigarIterator::new(cigars.into_iter());
        let mut results = Vec::new();
        while let Some(Ok((elem, count))) = collated.next() {
            results.push((elem.chrom_id, elem.reference_position, count));
        }
        // Should group by chrom_id and reference_position
        assert_eq!(results.len(), 2);
        assert_eq!(results[0], (1, 100, 2));
        assert_eq!(results[1], (2, 100, 2));
    }
    #[test]
    fn test_collated_augmented_cigar_iterator_multiple_same_position() {
        let cigars = vec![
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
            std::io::Result::Ok(("1M".to_string(), 1, 100)),
        ];
        let mut collated = CollatedAugmentedCigarIterator::new(cigars.into_iter());
        let mut results = Vec::new();
        while let Some(Ok((elem, count))) = collated.next() {
            results.push((elem, count));
        }
        // All three should be collated at position 100
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].0.reference_position, 100);
        assert_eq!(results[0].1, 3);
    }
}
