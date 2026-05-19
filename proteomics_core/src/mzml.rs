//! mzML spectrum parser.
//!
//! Streaming XML parse via quick-xml. Decodes base64 + optional zlib-compressed
//! binary arrays. Handles both 32-bit and 64-bit float encodings.
//!
//! HUPO-PSI CV terms used:
//! - MS:1000511 → ms level
//! - MS:1000016 → scan start time (minutes)
//! - MS:1000744 → selected ion m/z
//! - MS:1000041 → charge state
//! - MS:1000514 → m/z array
//! - MS:1000515 → intensity array
//! - MS:1000521 → 32-bit float
//! - MS:1000523 → 64-bit float
//! - MS:1000574 → zlib compression
//! - MS:1000576 → no compression

use std::io::Read;

use anyhow::{Context, Result};
use base64::{engine::general_purpose::STANDARD as B64, Engine};
use flate2::read::ZlibDecoder;
use quick_xml::events::Event;
use quick_xml::Reader;

use crate::types::Spectrum;

/// Parse all spectra from an mzML byte slice.
///
/// Returns MS1 and MS2 spectra in scan order. Large files are handled
/// efficiently because quick-xml is zero-copy for attribute bytes.
pub fn parse_mzml(data: &[u8]) -> Result<Vec<Spectrum>> {
    let mut reader = Reader::from_reader(data);
    reader.config_mut().trim_text(true);

    let mut spectra: Vec<Spectrum> = Vec::new();
    let mut current: Option<Spectrum> = None;

    // State for the current <binaryDataArray>.
    let mut in_bda = false;
    let mut bda_is_mz = false;
    let mut bda_is_intensity = false;
    let mut bda_zlib = false;
    let mut bda_32bit = true; // default: 32-bit float
    let mut bda_binary: Vec<u8> = Vec::new();
    let mut in_binary = false;

    let mut buf = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => match e.name().as_ref() {
                b"spectrum" => {
                    current = Some(Spectrum::default());
                }
                b"binaryDataArray" => {
                    in_bda = true;
                    bda_is_mz = false;
                    bda_is_intensity = false;
                    bda_zlib = false;
                    bda_32bit = true;
                    bda_binary.clear();
                }
                b"binary" => {
                    if in_bda {
                        in_binary = true;
                        bda_binary.clear();
                    }
                }
                _ => {}
            },

            Ok(Event::Empty(ref e)) => {
                // <cvParam .../> — most metadata lives here.
                if let Some(spec) = current.as_mut() {
                    let accession = attr_str(e, b"accession").unwrap_or_default();
                    let value = attr_str(e, b"value").unwrap_or_default();
                    match accession.as_str() {
                        "MS:1000511" => {
                            spec.ms_level = value.parse().unwrap_or(0);
                        }
                        "MS:1000016" => {
                            // Retention time in minutes → convert to seconds.
                            let rt_min: f32 = value.parse().unwrap_or(0.0);
                            spec.rt = rt_min * 60.0;
                        }
                        "MS:1000744" => {
                            spec.precursor_mz = value.parse().unwrap_or(0.0);
                        }
                        "MS:1000041" => {
                            spec.precursor_z = value.parse().unwrap_or(0);
                        }
                        _ => {}
                    }
                }
                if in_bda {
                    let accession = attr_str(e, b"accession").unwrap_or_default();
                    match accession.as_str() {
                        "MS:1000514" => bda_is_mz = true,
                        "MS:1000515" => bda_is_intensity = true,
                        "MS:1000574" => bda_zlib = true,
                        "MS:1000521" => bda_32bit = true,
                        "MS:1000523" => bda_32bit = false,
                        _ => {}
                    }
                }
            }

            Ok(Event::Text(ref e)) => {
                if in_binary && in_bda {
                    bda_binary.extend_from_slice(e.as_ref());
                }
            }

            Ok(Event::End(ref e)) => match e.name().as_ref() {
                b"binary" => {
                    in_binary = false;
                }
                b"binaryDataArray" => {
                    if let Some(spec) = current.as_mut() {
                        if bda_is_mz || bda_is_intensity {
                            match decode_binary(&bda_binary, bda_zlib, bda_32bit) {
                                Ok(vals) => {
                                    if bda_is_mz {
                                        spec.mz = vals;
                                    } else {
                                        spec.intensity = vals;
                                    }
                                }
                                Err(e) => {
                                    log::warn!("binary decode failed for scan {}: {e}", spec.scan);
                                }
                            }
                        }
                    }
                    in_bda = false;
                }
                b"spectrum" => {
                    if let Some(mut spec) = current.take() {
                        // Assign scan from index if not set by cvParam.
                        if spec.scan == 0 {
                            spec.scan = spectra.len() as u32 + 1;
                        }
                        spectra.push(spec);
                    }
                }
                _ => {}
            },

            Ok(Event::Eof) => break,
            Err(e) => return Err(anyhow::anyhow!("mzML parse error: {e}")),
            _ => {}
        }
        buf.clear();
    }

    Ok(spectra)
}

/// Decode a base64-encoded, optionally zlib-compressed binary array.
fn decode_binary(raw: &[u8], zlib: bool, is_32bit: bool) -> Result<Vec<f32>> {
    // Strip whitespace that can appear inside <binary>.
    let clean: Vec<u8> = raw
        .iter()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect();
    if clean.is_empty() {
        return Ok(vec![]);
    }

    let decoded = B64.decode(&clean).context("base64 decode")?;

    let bytes: Vec<u8> = if zlib {
        let mut dec = ZlibDecoder::new(&decoded[..]);
        let mut out = Vec::new();
        dec.read_to_end(&mut out).context("zlib decompress")?;
        out
    } else {
        decoded
    };

    if is_32bit {
        Ok(bytes
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect())
    } else {
        Ok(bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]) as f32)
            .collect())
    }
}

fn attr_str(e: &quick_xml::events::BytesStart, name: &[u8]) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.as_ref() == name)
        .and_then(|a| String::from_utf8(a.value.as_ref().to_vec()).ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_MZML: &[u8] = br#"<?xml version="1.0"?>
<mzML>
  <run>
    <spectrumList count="2">
      <spectrum index="0" id="scan=1" defaultArrayLength="3">
        <cvParam accession="MS:1000511" value="1"/>
        <cvParam accession="MS:1000016" value="0.5"/>
        <binaryDataArrayList count="2">
          <binaryDataArray>
            <cvParam accession="MS:1000514"/>
            <cvParam accession="MS:1000576"/>
            <cvParam accession="MS:1000521"/>
            <binary>AAAAQAAAAEAAAABI</binary>
          </binaryDataArray>
          <binaryDataArray>
            <cvParam accession="MS:1000515"/>
            <cvParam accession="MS:1000576"/>
            <cvParam accession="MS:1000521"/>
            <binary>AAAAQAAAAEAAAABi</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
      <spectrum index="1" id="scan=2" defaultArrayLength="0">
        <cvParam accession="MS:1000511" value="2"/>
        <cvParam accession="MS:1000016" value="1.0"/>
        <precursorList count="1">
          <precursor>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam accession="MS:1000744" value="450.23"/>
                <cvParam accession="MS:1000041" value="2"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
        </precursorList>
        <binaryDataArrayList count="2">
          <binaryDataArray>
            <cvParam accession="MS:1000514"/>
            <cvParam accession="MS:1000576"/>
            <cvParam accession="MS:1000521"/>
            <binary></binary>
          </binaryDataArray>
          <binaryDataArray>
            <cvParam accession="MS:1000515"/>
            <cvParam accession="MS:1000576"/>
            <cvParam accession="MS:1000521"/>
            <binary></binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"#;

    #[test]
    fn parse_ms_levels() {
        let spectra = parse_mzml(MINIMAL_MZML).expect("parse");
        assert_eq!(spectra.len(), 2);
        assert_eq!(spectra[0].ms_level, 1);
        assert_eq!(spectra[1].ms_level, 2);
    }

    #[test]
    fn parse_retention_time() {
        let spectra = parse_mzml(MINIMAL_MZML).expect("parse");
        // 0.5 minutes → 30 seconds
        assert!((spectra[0].rt - 30.0).abs() < 0.1);
        assert!((spectra[1].rt - 60.0).abs() < 0.1);
    }

    #[test]
    fn parse_precursor() {
        let spectra = parse_mzml(MINIMAL_MZML).expect("parse");
        assert!((spectra[1].precursor_mz - 450.23).abs() < 0.01);
        assert_eq!(spectra[1].precursor_z, 2);
    }
}
