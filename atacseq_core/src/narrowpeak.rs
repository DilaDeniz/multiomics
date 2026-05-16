use std::path::Path;

use anyhow::{Context, Result};
use memmap2::{Advice, Mmap, MmapOptions};

use biomics_core::parse::{parse_f64, parse_u64, ByteLines, TabFields};

use crate::types::AtacPeak;

/// Parse an ENCODE narrowPeak file and return a `Vec<AtacPeak>`.
///
/// The file is memory-mapped with `MADV_SEQUENTIAL` for efficient sequential
/// reads. Lines beginning with `#`, `track`, or `browser` are skipped.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or memory-mapped. Individual
/// malformed lines are logged as warnings and skipped rather than aborting.
pub fn parse_narrowpeak(path: &Path) -> Result<Vec<AtacPeak>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open narrowPeak file: {}", path.display()))?;

    // Safety: the file is opened read-only; we do not mutate the mapping.
    let mmap: Mmap = unsafe {
        MmapOptions::new()
            .map(&file)
            .with_context(|| format!("cannot mmap narrowPeak file: {}", path.display()))?
    };

    // Advise the kernel that we will read sequentially so it can pre-fetch.
    if let Err(e) = mmap.advise(Advice::Sequential) {
        log::debug!("madvise(SEQUENTIAL) failed (non-fatal): {}", e);
    }

    let data: &[u8] = &mmap;

    if data.is_empty() {
        return Ok(Vec::new());
    }

    // Heuristic: one line ≈ 80 bytes in a typical narrowPeak file.
    let mut peaks = Vec::with_capacity(data.len() / 80);

    for (line_no, line) in ByteLines::new(data).enumerate() {
        // Skip blank lines and comment / track-definition headers.
        if line.is_empty()
            || line.starts_with(b"#")
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }

        match parse_line(line) {
            Ok(peak) => peaks.push(peak),
            Err(e) => {
                log::warn!("narrowPeak line {}: {} (skipping)", line_no + 1, e);
            }
        }
    }

    Ok(peaks)
}

/// Parse a single narrowPeak data line (BED6+4 = 10 columns).
///
/// Returns an error when fewer than 10 tab-delimited columns are present or
/// any mandatory numeric field cannot be parsed.
#[inline]
fn parse_line(line: &[u8]) -> Result<AtacPeak> {
    let mut fields = TabFields::new(line);

    macro_rules! next_field {
        ($label:literal) => {{
            fields
                .next()
                .with_context(|| format!("missing field: {}", $label))?
        }};
    }

    // col 1 — chrom
    let chrom_bytes = next_field!("chrom");
    let chrom = std::str::from_utf8(chrom_bytes)
        .context("chrom is not valid UTF-8")?
        .to_owned();

    // col 2 — start (0-based)
    let start_bytes = next_field!("start");
    let start = parse_u64(start_bytes)
        .with_context(|| format!("invalid start: {}", lossy(start_bytes)))?;

    // col 3 — end
    let end_bytes = next_field!("end");
    let end =
        parse_u64(end_bytes).with_context(|| format!("invalid end: {}", lossy(end_bytes)))?;

    // col 4 — name
    let name_bytes = next_field!("name");
    let name = std::str::from_utf8(name_bytes)
        .context("name is not valid UTF-8")?
        .to_owned();

    // col 5 — score (integer 0–1000, stored as f64)
    let score_bytes = next_field!("score");
    let score = parse_f64(score_bytes)
        .with_context(|| format!("invalid score: {}", lossy(score_bytes)))?;

    // col 6 — strand: '+', '-', or '.' (dot → None)
    let strand_bytes = next_field!("strand");
    let strand = parse_strand(strand_bytes)
        .with_context(|| format!("invalid strand: {}", lossy(strand_bytes)))?;

    // col 7 — signalValue
    let sig_bytes = next_field!("signalValue");
    let signal_value = parse_f64(sig_bytes)
        .with_context(|| format!("invalid signalValue: {}", lossy(sig_bytes)))?;

    // col 8 — pValue (-log10; -1 if not computed)
    let pv_bytes = next_field!("pValue");
    let p_value_log10 = parse_f64(pv_bytes)
        .with_context(|| format!("invalid pValue: {}", lossy(pv_bytes)))?;

    // col 9 — qValue (-log10; -1 if not computed)
    let qv_bytes = next_field!("qValue");
    let q_value_log10 = parse_f64(qv_bytes)
        .with_context(|| format!("invalid qValue: {}", lossy(qv_bytes)))?;

    // col 10 — peak (0-based offset to summit; -1 if not determined)
    let pk_bytes = next_field!("peak");
    let peak_offset = parse_peak_offset(pk_bytes)
        .with_context(|| format!("invalid peak: {}", lossy(pk_bytes)))?;

    Ok(AtacPeak {
        chrom,
        start,
        end,
        name,
        score,
        strand,
        signal_value,
        p_value_log10,
        q_value_log10,
        peak_offset,
    })
}

/// Decode the strand column: `+` → `Some('+')`, `-` → `Some('-')`, `.` → `None`.
#[inline]
fn parse_strand(bytes: &[u8]) -> Result<Option<char>> {
    match bytes {
        b"+" => Ok(Some('+')),
        b"-" => Ok(Some('-')),
        b"." | b"" => Ok(None),
        other => anyhow::bail!("unexpected strand value: {}", lossy(other)),
    }
}

/// Decode the peak-offset column. The spec uses `-1` to signal "not computed".
/// We support both `-1` (i64 sentinel) and any non-negative offset.
#[inline]
fn parse_peak_offset(bytes: &[u8]) -> Result<i64> {
    // Fast path: signed integer; the only negative value in the spec is -1.
    if bytes == b"-1" {
        return Ok(-1);
    }
    let n = parse_u64(bytes).with_context(|| format!("expected integer, got {}", lossy(bytes)))?;
    Ok(n as i64)
}

/// Lossily convert bytes to a `String` for error messages (no heap alloc path
/// for valid ASCII, but this is only called on the error path).
#[inline]
fn lossy(b: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(b)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Write bytes to a temp file and return its path (as a `NamedTempFile` so it
    /// is cleaned up automatically).
    fn temp_narrowpeak(content: &[u8]) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(content).expect("write");
        f.flush().expect("flush");
        f
    }

    #[test]
    fn test_empty_file() {
        let f = temp_narrowpeak(b"");
        let peaks = parse_narrowpeak(f.path()).expect("parse");
        assert!(peaks.is_empty(), "expected no peaks from empty file");
    }

    #[test]
    fn test_valid_ten_column_line() {
        let line = b"chr1\t1000\t2000\tpeak_001\t500\t.\t8.5\t6.2\t4.1\t450\n";
        let f = temp_narrowpeak(line);
        let peaks = parse_narrowpeak(f.path()).expect("parse");
        assert_eq!(peaks.len(), 1);

        let p = &peaks[0];
        assert_eq!(p.chrom, "chr1");
        assert_eq!(p.start, 1000);
        assert_eq!(p.end, 2000);
        assert_eq!(p.name, "peak_001");
        assert!((p.score - 500.0).abs() < f64::EPSILON);
        assert_eq!(p.strand, None);
        assert!((p.signal_value - 8.5).abs() < 1e-10);
        assert!((p.p_value_log10 - 6.2).abs() < 1e-10);
        assert!((p.q_value_log10 - 4.1).abs() < 1e-10);
        assert_eq!(p.peak_offset, 450);
    }

    #[test]
    fn test_truncated_six_column_line_is_skipped() {
        // Only 6 columns — parse_line should fail and the line should be skipped.
        let line = b"chr2\t5000\t6000\tpeak_002\t300\t+\n";
        let f = temp_narrowpeak(line);
        // parse_narrowpeak logs a warning but must not return an Err itself.
        let peaks = parse_narrowpeak(f.path()).expect("parse");
        assert!(
            peaks.is_empty(),
            "truncated line should be skipped, got {} peaks",
            peaks.len()
        );
    }

    #[test]
    fn test_comment_and_track_lines_skipped() {
        let content = b"# comment\ntrack name=test\nbrowser position chr1:1-100\nchr3\t100\t200\tp1\t900\t+\t12.0\t8.0\t5.0\t50\n";
        let f = temp_narrowpeak(content);
        let peaks = parse_narrowpeak(f.path()).expect("parse");
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].chrom, "chr3");
    }

    #[test]
    fn test_peak_offset_minus_one_sentinel() {
        let line = b"chrX\t0\t100\tp\t0\t.\t1.0\t-1\t-1\t-1\n";
        let f = temp_narrowpeak(line);
        let peaks = parse_narrowpeak(f.path()).expect("parse");
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].peak_offset, -1);
        assert!((peaks[0].p_value_log10 - (-1.0)).abs() < f64::EPSILON);
    }
}
