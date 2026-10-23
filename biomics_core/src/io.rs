//! Input helpers shared by the text-format parsers (VCF, BED, TSV, narrowPeak).
//!
//! Plain files are memory-mapped so parsing stays zero-copy. Files that start
//! with the gzip magic bytes (`1f 8b`) — plain gzip as well as the BGZF blocks
//! written by `bgzip`/`bcftools`/`tabix` — are decompressed into memory first.
//! BGZF is a series of concatenated gzip members, so a multi-member decoder is
//! required; a single-member decoder would stop after the first 64 KiB block.

use std::io::Read;
use std::ops::Deref;
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use memmap2::Mmap;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// The bytes of an input file, either memory-mapped or decompressed.
pub enum InputBytes {
    /// Uncompressed file mapped read-only into memory.
    Mapped(Mmap),
    /// Decompressed (or empty) file contents held on the heap.
    Owned(Vec<u8>),
}

impl Deref for InputBytes {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        match self {
            InputBytes::Mapped(m) => m,
            InputBytes::Owned(v) => v,
        }
    }
}

impl AsRef<[u8]> for InputBytes {
    fn as_ref(&self) -> &[u8] {
        self
    }
}

/// Returns true when `data` begins with the gzip magic bytes.
#[inline]
pub fn is_gzip(data: &[u8]) -> bool {
    data.starts_with(&GZIP_MAGIC)
}

/// Open `path` for parsing, transparently decompressing gzip/BGZF input.
pub fn open_input(path: &Path) -> Result<InputBytes> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open '{}'", path.display()))?;
    let len = file
        .metadata()
        .with_context(|| format!("Cannot stat '{}'", path.display()))?
        .len();
    if len == 0 {
        return Ok(InputBytes::Owned(Vec::new()));
    }

    // SAFETY: the mapping is read-only and the file is not modified while the
    // process holds it.
    let mmap =
        unsafe { Mmap::map(&file) }.with_context(|| format!("Cannot mmap '{}'", path.display()))?;

    if !is_gzip(&mmap) {
        let _ = mmap.advise(memmap2::Advice::Sequential);
        return Ok(InputBytes::Mapped(mmap));
    }

    let mut out = Vec::with_capacity(mmap.len().saturating_mul(4));
    MultiGzDecoder::new(&mmap[..])
        .read_to_end(&mut out)
        .with_context(|| format!("Cannot decompress gzip input '{}'", path.display()))?;
    Ok(InputBytes::Owned(out))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        enc.finish().unwrap()
    }

    fn temp_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("biomics_io_{}_{}", std::process::id(), tag))
    }

    #[test]
    fn plain_file_is_mapped() {
        let path = temp_path("plain.txt");
        std::fs::write(&path, b"chr1\t1\n").unwrap();
        let bytes = open_input(&path).unwrap();
        assert!(matches!(bytes, InputBytes::Mapped(_)));
        assert_eq!(&*bytes, b"chr1\t1\n");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn multi_member_gzip_is_fully_decompressed() {
        // BGZF files are concatenated gzip members; all must be read.
        let mut data = gzip(b"line1\n");
        data.extend(gzip(b"line2\n"));
        let path = temp_path("multi.gz");
        std::fs::write(&path, &data).unwrap();
        let bytes = open_input(&path).unwrap();
        assert_eq!(&*bytes, b"line1\nline2\n");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn empty_file_yields_no_bytes() {
        let path = temp_path("empty");
        std::fs::write(&path, b"").unwrap();
        assert!(open_input(&path).unwrap().is_empty());
        std::fs::remove_file(path).unwrap();
    }
}
