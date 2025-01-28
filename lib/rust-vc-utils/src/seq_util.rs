pub fn rev_comp_in_place(dna: &mut [u8]) {
    let comp_base = |x: &mut u8| {
        *x = match *x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => *x,
        };
    };

    let len = dna.len();
    let halflen = len - len / 2;
    for i in 0..halflen {
        comp_base(&mut dna[i]);
        let rev_i = len - 1 - i;
        if i != rev_i {
            comp_base(&mut dna[rev_i]);
            dna.swap(i, rev_i);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rev_comp_in_place() {
        let mut input = b"NNATGCG".to_vec();
        let expected_output = b"CGCATNN".to_vec();
        rev_comp_in_place(&mut input);
        assert_eq!(input, expected_output);
    }
}
