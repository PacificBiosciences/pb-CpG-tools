use std::collections::BTreeMap;

/// Very simple ring buffer implementation
#[derive(Clone)]
pub struct RingBuffer<T> {
    max_size: usize,
    head_index: usize,
    data: Vec<T>,
}

impl<T> RingBuffer<T> {
    pub fn new(max_size: usize) -> Self {
        Self {
            max_size,
            head_index: 0,
            data: Vec::new(),
        }
    }

    /// Translate user index to vector index
    fn get_vec_index(&self, index: usize) -> usize {
        assert!(index < self.max_size);
        (index + self.head_index) % self.max_size
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn get_item(&self, index: usize) -> &T {
        assert!(index < self.data.len());
        &self.data[self.get_vec_index(index)]
    }

    pub fn push_item(&mut self, item: T) {
        if self.data.len() < self.max_size {
            self.data.push(item);
        } else {
            let vec_index = self.get_vec_index(0);
            self.data[vec_index] = item;
            self.head_index = (self.head_index + 1) % self.max_size
        }
    }
}

/// Track the total content within a fixed window_size, assuming observations at
/// monotonically increasing positions
///
pub struct SparseWindowSum {
    window_size: usize,
    sum: u32,
    pos_map: BTreeMap<i64, u32>,
}

impl SparseWindowSum {
    pub fn new(window_size: usize) -> Self {
        assert!(window_size > 1);
        Self {
            window_size,
            sum: 0,
            pos_map: BTreeMap::new(),
        }
    }

    pub fn sum(&self) -> u32 {
        self.sum
    }

    pub fn clear(&mut self) {
        self.sum = 0;
        self.pos_map.clear();
    }

    /// Index must be equal or higher than any previous value:
    pub fn push(&mut self, pos: i64, count: u32) {
        let mut clear = false;
        if let Some((last_pos, _last_count)) = self.pos_map.last_key_value() {
            assert!(pos > *last_pos);
            if (pos - *last_pos) >= self.window_size as i64 {
                clear = true;
            }
        }

        if clear {
            self.clear();
        } else {
            // Iterate through current members and trim every key below min_pos:
            let min_pos = 1 + pos - self.window_size as i64;
            let mut remove_pos = Vec::new();
            for (map_pos, map_count) in self.pos_map.iter() {
                if *map_pos >= min_pos {
                    break;
                }
                assert!(*map_count <= self.sum);
                self.sum -= *map_count;
                remove_pos.push(*map_pos);
            }

            for map_pos in remove_pos {
                self.pos_map.remove(&map_pos);
            }
        }

        self.pos_map.insert(pos, count);
        self.sum += count;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ring_buffer() {
        let mut rb = RingBuffer::new(3);

        assert_eq!(0, rb.len());

        rb.push_item(22);
        rb.push_item(8);

        assert_eq!(2, rb.len());
        assert_eq!(&8, rb.get_item(1));

        rb.push_item(6);
        rb.push_item(1);

        assert_eq!(3, rb.len());
        assert_eq!(&6, rb.get_item(1));
    }

    #[test]
    fn test_sparse_window_sum() {
        let mut sws = SparseWindowSum::new(3);

        assert_eq!(0, sws.sum());

        sws.push(100, 2);
        assert_eq!(2, sws.sum());
        sws.push(101, 2);
        assert_eq!(4, sws.sum());
        sws.push(102, 2);
        assert_eq!(6, sws.sum());
        sws.push(103, 2);
        assert_eq!(6, sws.sum());

        sws.push(200, 2);
        assert_eq!(2, sws.sum());
    }
}
