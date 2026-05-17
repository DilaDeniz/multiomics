use aho_corasick::AhoCorasick;
use memchr::memchr;

/// Zero-allocation line iterator over a byte slice.
///
/// Yields `&[u8]` slices, stripping the trailing `\r` on Windows line endings.
/// Avoids the `String` allocation that `BufReader::lines()` performs on every line.
pub struct ByteLines<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> ByteLines<'a> {
    #[inline]
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }
}

impl<'a> Iterator for ByteLines<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.data.len() {
            return None;
        }
        let start = self.pos;
        match memchr(b'\n', &self.data[start..]) {
            Some(n) => {
                let end = start + n;
                self.pos = end + 1;
                let line = &self.data[start..end];
                Some(if line.last() == Some(&b'\r') {
                    &line[..line.len() - 1]
                } else {
                    line
                })
            }
            None => {
                self.pos = self.data.len();
                Some(&self.data[start..])
            }
        }
    }
}

/// Zero-allocation tab-field iterator over a byte slice.
///
/// Yields each `\t`-delimited field as a `&[u8]`. The final field is yielded
/// even when it is not followed by a tab.
pub struct TabFields<'a> {
    data: &'a [u8],
    pos: usize,
    done: bool,
}

impl<'a> TabFields<'a> {
    #[inline]
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            pos: 0,
            done: false,
        }
    }
}

impl<'a> Iterator for TabFields<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let start = self.pos;
        if start > self.data.len() {
            return None;
        }
        match memchr(b'\t', &self.data[start..]) {
            Some(n) => {
                let end = start + n;
                self.pos = end + 1;
                Some(&self.data[start..end])
            }
            None => {
                self.done = true;
                Some(&self.data[start..])
            }
        }
    }
}

/// Trim leading and trailing ASCII whitespace from a byte slice.
#[inline]
pub fn trim_bytes(b: &[u8]) -> &[u8] {
    let start = b
        .iter()
        .position(|&c| c != b' ' && c != b'\t' && c != b'\r' && c != b'\n')
        .unwrap_or(b.len());
    let end = b
        .iter()
        .rposition(|&c| c != b' ' && c != b'\t' && c != b'\r' && c != b'\n')
        .map(|i| i + 1)
        .unwrap_or(0);
    if start < end {
        &b[start..end]
    } else {
        b""
    }
}

/// Parse a decimal `u64` from bytes. Returns `None` if any byte is non-ASCII-digit.
///
/// ~5x faster than `str::parse::<u64>()` for typical genomic coordinate strings.
#[inline]
pub fn parse_u64(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() {
        return None;
    }
    let mut n = 0u64;
    for &b in bytes {
        let d = b.wrapping_sub(b'0');
        if d > 9 {
            return None;
        }
        n = n.wrapping_mul(10).wrapping_add(d as u64);
    }
    Some(n)
}

/// Parse a `f64` from bytes via `fast-float`.
///
/// Roughly 5–10× faster than `str::parse::<f64>()` on typical floating-point
/// strings like TPM values or quality scores.
#[inline]
pub fn parse_f64(bytes: &[u8]) -> Option<f64> {
    fast_float::parse(bytes).ok()
}

/// Parse a `f32` from bytes via `fast-float`.
#[inline]
pub fn parse_f32(bytes: &[u8]) -> Option<f32> {
    fast_float::parse(bytes).ok()
}

/// Find the n-th `|`-delimited field (0-indexed) in a byte slice.
#[inline]
pub fn nth_pipe_field(data: &[u8], n: usize) -> Option<&[u8]> {
    let mut count = 0usize;
    let mut start = 0usize;
    for (i, &b) in data.iter().enumerate() {
        if b == b'|' {
            if count == n {
                return Some(&data[start..i]);
            }
            count += 1;
            start = i + 1;
        }
    }
    if count == n {
        Some(&data[start..])
    } else {
        None
    }
}

/// Search `info` (a `KEY=val;KEY2=val2;…` byte slice) for `key` and return
/// the value slice. Avoids any heap allocation.
#[inline]
pub fn info_value_bytes<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    let klen = key.len();
    let mut pos = 0usize;
    while pos < info.len() {
        let field_end = memchr(b';', &info[pos..])
            .map(|n| pos + n)
            .unwrap_or(info.len());
        let field = &info[pos..field_end];
        if field.len() > klen && &field[..klen] == key && field[klen] == b'=' {
            return Some(&field[klen + 1..]);
        }
        pos = field_end + 1;
    }
    None
}

/// Single-pass multi-key extractor for VCF INFO fields.
///
/// Builds an aho-corasick automaton over the supplied `key=` byte patterns and
/// scans the INFO string exactly once, returning one `Option<&[u8]>` per key.
/// Correct field-boundary semantics: a key is only accepted when it appears at
/// position 0 or immediately after a `;` separator.
///
/// # Performance
/// For k independent keys, the naïve approach calls `info_value_bytes` k times
/// → k passes over the INFO string. `InfoMultiParser` reduces this to one pass
/// regardless of k, giving a k× speedup for INFO-heavy parsers (e.g. CNV VCF).
///
/// # Example
/// ```
/// use biomics_core::parse::InfoMultiParser;
/// static PARSER: std::sync::LazyLock<InfoMultiParser> =
///     std::sync::LazyLock::new(|| InfoMultiParser::new(&["SVTYPE", "CN", "END"]));
/// let vals = PARSER.extract(b"SVTYPE=DEL;CN=1;END=5000");
/// assert_eq!(vals[0], Some(b"DEL".as_ref()));
/// assert_eq!(vals[1], Some(b"1".as_ref()));
/// assert_eq!(vals[2], Some(b"5000".as_ref()));
/// ```
pub struct InfoMultiParser {
    ac: AhoCorasick,
    n_keys: usize,
}

impl InfoMultiParser {
    /// Build the automaton. `keys` are plain key names (without the `=`
    /// suffix); the parser appends `=` internally.
    pub fn new(keys: &[&str]) -> Self {
        let patterns: Vec<Vec<u8>> = keys.iter().map(|k| format!("{k}=").into_bytes()).collect();
        let ac = AhoCorasick::builder()
            .ascii_case_insensitive(false)
            .build(&patterns)
            .expect("InfoMultiParser: invalid patterns");
        Self {
            ac,
            n_keys: keys.len(),
        }
    }

    /// Extract all key values from one INFO byte slice in a single scan.
    ///
    /// Returns a `Vec<Option<&[u8]>>` parallel to the `keys` slice passed to
    /// [`InfoMultiParser::new`]. Fields are accepted only at field-start
    /// positions (byte 0 or immediately after `;`).
    pub fn extract<'a>(&self, info: &'a [u8]) -> Vec<Option<&'a [u8]>> {
        let mut out = vec![None; self.n_keys];
        let mut remaining = self.n_keys;

        for mat in self.ac.find_iter(info) {
            let start = mat.start();
            // Enforce field-start: must be at byte 0 or after ';'
            if start != 0 && info.get(start - 1) != Some(&b';') {
                continue;
            }
            let key_idx = mat.pattern().as_usize();
            if out[key_idx].is_some() {
                continue; // already captured — take first occurrence
            }
            let val_start = mat.end(); // byte right after 'key='
            let val_end = memchr(b';', &info[val_start..])
                .map(|n| val_start + n)
                .unwrap_or(info.len());
            out[key_idx] = Some(&info[val_start..val_end]);
            remaining -= 1;
            if remaining == 0 {
                break;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_byte_lines() {
        let data = b"line1\nline2\r\nline3";
        let lines: Vec<_> = ByteLines::new(data).collect();
        assert_eq!(lines, [b"line1".as_ref(), b"line2", b"line3"]);
    }

    #[test]
    fn test_tab_fields() {
        let data = b"chr1\t100\t200\t85.0";
        let fields: Vec<_> = TabFields::new(data).collect();
        assert_eq!(fields, [b"chr1".as_ref(), b"100", b"200", b"85.0"]);
    }

    #[test]
    fn test_parse_u64() {
        assert_eq!(parse_u64(b"12345678"), Some(12345678));
        assert_eq!(parse_u64(b"0"), Some(0));
        assert_eq!(parse_u64(b"abc"), None);
        assert_eq!(parse_u64(b""), None);
    }

    #[test]
    fn test_info_value_bytes() {
        let info = b"AF=0.42;DP=30;GENE=TP53";
        assert_eq!(info_value_bytes(info, b"AF"), Some(b"0.42".as_ref()));
        assert_eq!(info_value_bytes(info, b"GENE"), Some(b"TP53".as_ref()));
        assert_eq!(info_value_bytes(info, b"MISSING"), None);
    }
}
