use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;

const SVG_SIZE: f64 = 600.0;
const CX: f64 = SVG_SIZE / 2.0;
const CY: f64 = SVG_SIZE / 2.0;

const R_OUTER: f64 = 240.0;
const R_MIDDLE: f64 = 195.0;
const R_INNER: f64 = 150.0;
const TRACK_W: f64 = 38.0;

const CHROM_SIZES_MB: &[(&str, f64)] = &[
    ("chr1", 248.9),
    ("chr2", 242.2),
    ("chr3", 198.3),
    ("chr4", 190.2),
    ("chr5", 181.5),
    ("chr6", 170.8),
    ("chr7", 159.3),
    ("chr8", 145.1),
    ("chr9", 138.4),
    ("chr10", 133.8),
    ("chr11", 135.1),
    ("chr12", 133.3),
    ("chr13", 114.4),
    ("chr14", 107.0),
    ("chr15", 101.9),
    ("chr16", 90.3),
    ("chr17", 83.3),
    ("chr18", 80.4),
    ("chr19", 58.6),
    ("chr20", 64.4),
    ("chr21", 46.7),
    ("chr22", 50.8),
    ("chrX", 156.0),
    ("chrY", 57.2),
];

pub fn generate_circos_svg(genomics: &GenomicsSummary, epigen: &EpigenomicsSummary) -> String {
    let total_genome_mb: f64 = CHROM_SIZES_MB.iter().map(|(_, s)| s).sum();

    let max_variants = genomics
        .per_chrom
        .values()
        .map(|d| d.total)
        .max()
        .unwrap_or(1)
        .max(1) as f64;

    let gap_deg = 1.5_f64;
    let total_gap = gap_deg * CHROM_SIZES_MB.len() as f64;
    let data_deg = 360.0 - total_gap;

    let mut s = String::with_capacity(32_768);

    let sz = SVG_SIZE as u32;
    let cx = CX as u32;

    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {sz} {sz}\" \
         width=\"{sz}\" height=\"{sz}\" style=\"background:#0d1117;font-family:'Segoe UI',sans-serif\">\n\
         <defs><filter id=\"glow\"><feGaussianBlur stdDeviation=\"2\" result=\"blur\"/>\
         <feMerge><feMergeNode in=\"blur\"/><feMergeNode in=\"SourceGraphic\"/>\
         </feMerge></filter></defs>\n\
         <text x=\"{cx}\" y=\"28\" text-anchor=\"middle\" fill=\"#58a6ff\" \
         font-size=\"14\" font-weight=\"600\">Genomic Overview</text>\n"
    ));

    s.push_str(
        "<rect x=\"10\" y=\"40\" width=\"12\" height=\"12\" fill=\"#3fb950\" rx=\"2\"/>\
         <text x=\"26\" y=\"51\" fill=\"#8b949e\" font-size=\"11\">Variant density</text>\n\
         <rect x=\"10\" y=\"58\" width=\"12\" height=\"12\" fill=\"#bc8cff\" rx=\"2\"/>\
         <text x=\"26\" y=\"69\" fill=\"#8b949e\" font-size=\"11\">Methylation</text>\n\
         <rect x=\"10\" y=\"76\" width=\"12\" height=\"12\" fill=\"#161b22\" stroke=\"#30363d\" rx=\"2\"/>\
         <text x=\"26\" y=\"87\" fill=\"#8b949e\" font-size=\"11\">Chromosomes</text>\n"
    );

    let mut angle_deg = -90.0_f64;

    for (chrom, size_mb) in CHROM_SIZES_MB {
        let arc_deg = size_mb / total_genome_mb * data_deg;
        let start_rad = angle_deg.to_radians();
        let end_rad = (angle_deg + arc_deg).to_radians();
        let mid_rad = (angle_deg + arc_deg / 2.0).to_radians();
        let large = if arc_deg > 180.0 { 1 } else { 0 };

        // Chromosome backbone ring
        let ri = R_INNER;
        let ro = R_INNER + 12.0;
        let x1 = CX + ri * start_rad.cos();
        let y1 = CY + ri * start_rad.sin();
        let x2 = CX + ro * start_rad.cos();
        let y2 = CY + ro * start_rad.sin();
        let x3 = CX + ro * end_rad.cos();
        let y3 = CY + ro * end_rad.sin();
        let x4 = CX + ri * end_rad.cos();
        let y4 = CY + ri * end_rad.sin();

        s.push_str(&format!(
            "<path d=\"M {x1:.2} {y1:.2} A {ri:.2} {ri:.2} 0 {large} 1 {x4:.2} {y4:.2} \
             L {x3:.2} {y3:.2} A {ro:.2} {ro:.2} 0 {large} 0 {x2:.2} {y2:.2} Z\" \
             fill=\"#21262d\" stroke=\"#30363d\" stroke-width=\"0.5\"/>\n"
        ));

        if arc_deg > 6.0 {
            let label_r = R_INNER - 18.0;
            let lx = CX + label_r * mid_rad.cos();
            let ly = CY + label_r * mid_rad.sin();
            let label = chrom.strip_prefix("chr").unwrap_or(chrom);
            let rot = (angle_deg + arc_deg / 2.0 + 90.0) % 360.0;
            s.push_str(&format!(
                "<text x=\"{lx:.1}\" y=\"{ly:.1}\" fill=\"#8b949e\" font-size=\"8\" \
                 text-anchor=\"middle\" dominant-baseline=\"middle\" \
                 transform=\"rotate({rot:.1},{lx:.1},{ly:.1})\">{label}</text>\n"
            ));
        }

        // Methylation ring (middle)
        let meth_pct = epigen
            .per_chrom
            .get(*chrom)
            .map(|cm| cm.mean_methylation)
            .unwrap_or(75.0);
        let meth_intensity = (meth_pct / 100.0).clamp(0.0, 1.0);
        let meth_color = interpolate_color(0x58, 0xa6, 0xff, 0xbc, 0x8c, 0xff, meth_intensity);
        push_ring_segment(
            &mut s,
            R_MIDDLE,
            R_MIDDLE + TRACK_W,
            start_rad,
            end_rad,
            arc_deg,
            &meth_color,
        );

        // Variant density ring (outer)
        let variants = genomics
            .per_chrom
            .get(*chrom)
            .map(|d| d.total as f64)
            .unwrap_or(0.0);
        let density = (variants / max_variants).clamp(0.0, 1.0);
        let track_h = (density * TRACK_W).max(2.0);
        let var_color = interpolate_color(0x23, 0x48, 0x23, 0x3f, 0xb9, 0x50, density);
        push_ring_segment(
            &mut s,
            R_OUTER,
            R_OUTER + track_h,
            start_rad,
            end_rad,
            arc_deg,
            &var_color,
        );

        angle_deg += arc_deg + gap_deg;
    }

    let cy1 = (CY - 8.0) as u32;
    let cy2 = (CY + 10.0) as u32;
    let tv = format!("{:.1}M vars", genomics.total_variants as f64 / 1_000_000.0);
    let mv = epigen.global_methylation_pct;

    s.push_str(&format!(
        "<text x=\"{cx}\" y=\"{cy1}\" text-anchor=\"middle\" fill=\"#c9d1d9\" \
         font-size=\"13\" font-weight=\"600\">{tv}</text>\n\
         <text x=\"{cx}\" y=\"{cy2}\" text-anchor=\"middle\" fill=\"#8b949e\" \
         font-size=\"10\">{mv:.1}% meth</text>\n"
    ));

    s.push_str("</svg>");
    s
}

// cx/cy are always CX/CY constants — callers omit them to stay under the 7-arg limit.
fn push_ring_segment(
    out: &mut String,
    r_inner: f64,
    r_outer: f64,
    start: f64,
    end: f64,
    arc_deg: f64,
    fill: &str,
) {
    let large = if arc_deg > 180.0 { 1 } else { 0 };
    let x1 = CX + r_inner * start.cos();
    let y1 = CY + r_inner * start.sin();
    let x2 = CX + r_outer * start.cos();
    let y2 = CY + r_outer * start.sin();
    let x3 = CX + r_outer * end.cos();
    let y3 = CY + r_outer * end.sin();
    let x4 = CX + r_inner * end.cos();
    let y4 = CY + r_inner * end.sin();

    out.push_str(&format!(
        "<path d=\"M {x1:.2} {y1:.2} A {ri:.2} {ri:.2} 0 {large} 1 {x4:.2} {y4:.2} \
         L {x3:.2} {y3:.2} A {ro:.2} {ro:.2} 0 {large} 0 {x2:.2} {y2:.2} Z\" \
         fill=\"{fill}\" opacity=\"0.85\"/>\n",
        ri = r_inner,
        ro = r_outer,
    ));
}

fn interpolate_color(r1: u8, g1: u8, b1: u8, r2: u8, g2: u8, b2: u8, t: f64) -> String {
    let r = (r1 as f64 + (r2 as f64 - r1 as f64) * t) as u8;
    let g = (g1 as f64 + (g2 as f64 - g1 as f64) * t) as u8;
    let b = (b1 as f64 + (b2 as f64 - b1 as f64) * t) as u8;
    format!("#{r:02X}{g:02X}{b:02X}")
}
