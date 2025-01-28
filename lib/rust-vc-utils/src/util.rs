/// Updatable/mergable mean value tracker
///
#[derive(Default, Clone)]
pub struct MeanTracker {
    total: f64,
    count: f64,
}

impl MeanTracker {
    pub fn mean(&self) -> f64 {
        if self.count > 0.0 {
            self.total / self.count
        } else {
            0.0
        }
    }

    pub fn insert(&mut self, x: f64) {
        self.total += x;
        self.count += 1.0;
    }

    pub fn merge(&mut self, other: &Self) {
        self.total += other.total;
        self.count += other.count;
    }
}

/// Deterministically downsample a vector while evenly distributing the removed positions
///
pub fn downsample_vector<T>(vec: Vec<T>, new_size: usize) -> Vec<T> {
    let size = vec.len();
    if size <= new_size {
        vec
    } else {
        let mut t = 0;
        vec.into_iter()
            .filter(|_| {
                t = (t % size) + new_size;
                t >= size
            })
            .collect()
    }
}

/// Divide a region into segments so that no segment is larger than segment_size
///
/// Return a series of zero-indexed half-closed (begin,end) intervals defining the segments
///
pub fn get_region_segments(size: u64, segment_size: u64) -> Vec<(u64, u64)> {
    let segment_count = 1 + ((size - 1) / segment_size);
    let segment_base_size = size / segment_count;
    let n_plus_one = size % segment_count;

    let mut intervals = Vec::new();
    let mut start = 0;
    for segment_index in 0..segment_count {
        let mut segment_size = segment_base_size;
        if segment_index < n_plus_one {
            segment_size += 1;
        }
        let end = std::cmp::min(start + segment_size, size);
        intervals.push((start, end));
        start = end;
    }
    intervals
}

/// Creates an iterator for array, that yields bin ranges of non-excluded regions, given a specified
/// exclusion function.
///
/// # Example
///
/// If array is length 10 and item 3 (zero-indexed) has is excluded, this object yields the ranges
/// (0,3) and (4,10) when iterated.
///
pub struct ArraySegmenter<'a, T, U>
where
    U: Fn(&T) -> bool,
{
    current: std::ops::Range<usize>,
    array: &'a [T],
    exclude_func: U,
}

impl<'a, T, U> ArraySegmenter<'a, T, U>
where
    U: Fn(&T) -> bool,
{
    ///
    /// # Arguments
    ///
    /// * `exclude_func` - Function that returns true for any array element that should be excluded
    ///
    pub fn new(array: &'a [T], exclude_func: U) -> Self {
        Self {
            current: 0..0,
            array,
            exclude_func,
        }
    }
}

impl<'a, T, U> Iterator for ArraySegmenter<'a, T, U>
where
    U: Fn(&T) -> bool,
{
    type Item = std::ops::Range<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        let size = self.array.len();
        if self.current.end >= size {
            return None;
        }

        // Update start to point to the first non-excluded bin not less than end
        for bin in self.current.end..=size {
            if bin < size && (self.exclude_func)(&self.array[bin]) {
                continue;
            }
            self.current.start = bin;
            break;
        }

        // Update end to point to the first excluded bin not less than start
        for bin in self.current.start..=size {
            if bin < size && !(self.exclude_func)(&self.array[bin]) {
                continue;
            }
            self.current.end = bin;
            break;
        }
        Some(self.current.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean_tracker() {
        let mut a = MeanTracker::default();
        let mut b = MeanTracker::default();

        a.insert(2.0);
        a.insert(4.0);

        approx::assert_ulps_eq!(a.mean(), 3.0, max_ulps = 4);

        b.insert(6.0);
        a.merge(&b);

        approx::assert_ulps_eq!(a.mean(), 4.0, max_ulps = 4);
    }

    #[test]
    fn test_downsample_vector() {
        let t = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

        for i in 0..15 {
            assert_eq!(downsample_vector(t.clone(), i).len(), std::cmp::min(i, 13));
        }
    }

    #[test]
    fn test_get_region_segments() {
        assert_eq!(get_region_segments(100, 200), vec![(0, 100)]);
        assert_eq!(
            get_region_segments(100, 49),
            vec![(0, 34), (34, 67), (67, 100)]
        );
    }

    #[test]
    fn test_array_segmenter() {
        let test_input = vec![0, 1, 2, -1, 4, 5, 6, 7, 8, 9];
        let ranges = ArraySegmenter::new(&test_input, |x| *x < 0).collect::<Vec<_>>();

        assert_eq!(ranges.len(), 2);
        assert_eq!(ranges[0], 0..3);
        assert_eq!(ranges[1], 4..10);
    }
}
