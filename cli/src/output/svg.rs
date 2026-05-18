/// Lightweight SVG builder — collects SVG elements as a String.
pub struct Svg {
    width: u32,
    height: u32,
    elements: Vec<String>,
}

impl Svg {
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            elements: Vec::new(),
        }
    }

    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str, opacity: f64) -> &mut Self {
        self.elements.push(format!(
            r#"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="{fill}" opacity="{opacity:.2}"/>"#
        ));
        self
    }

    pub fn line(
        &mut self,
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        color: &str,
        width: f64,
    ) -> &mut Self {
        self.elements.push(format!(
            r#"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="{color}" stroke-width="{width:.2}"/>"#
        ));
        self
    }

    pub fn circle(&mut self, cx: f64, cy: f64, r: f64, fill: &str, opacity: f64) -> &mut Self {
        self.elements.push(format!(
            r#"<circle cx="{cx:.2}" cy="{cy:.2}" r="{r:.2}" fill="{fill}" opacity="{opacity:.2}"/>"#
        ));
        self
    }

    pub fn text(&mut self, x: f64, y: f64, s: &str, size: u32, anchor: &str) -> &mut Self {
        let escaped = escape_svg(s);
        self.elements.push(format!(
            r#"<text x="{x:.2}" y="{y:.2}" font-size="{size}" text-anchor="{anchor}" font-family="Arial, sans-serif">{escaped}</text>"#
        ));
        self
    }

    /// Returns a complete `<svg>…</svg>` string.
    pub fn build(self) -> String {
        let body = self.elements.join("\n");
        format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">{body}</svg>"#,
            w = self.width,
            h = self.height,
        )
    }
}

// ── Layout constants ───────────────────────────────────────────────────────────

const MARGIN_LEFT: f64 = 60.0;
const MARGIN_TOP: f64 = 40.0;
const MARGIN_BOTTOM: f64 = 60.0;
const MARGIN_RIGHT: f64 = 20.0;

/// Bar chart: title at top, vertical bars with labels below, Y-axis with 5 ticks & gridlines.
pub fn bar_chart_svg(
    title: &str,
    labels: &[&str],
    values: &[f64],
    color: &str,
    width: u32,
    height: u32,
) -> String {
    let mut svg = Svg::new(width, height);
    // White background
    svg.rect(0.0, 0.0, width as f64, height as f64, "#ffffff", 1.0);

    let pw = width as f64 - MARGIN_LEFT - MARGIN_RIGHT;
    let ph = height as f64 - MARGIN_TOP - MARGIN_BOTTOM;
    let ox = MARGIN_LEFT;
    let oy = MARGIN_TOP;

    // Title
    svg.text(ox + pw / 2.0, oy - 10.0, title, 14, "middle");

    if values.is_empty() {
        return svg.build();
    }

    let max_val = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let max_val = if max_val <= 0.0 { 1.0 } else { max_val };

    // Y-axis gridlines & ticks (5 intervals)
    let n_ticks = 5u32;
    for i in 0..=n_ticks {
        let frac = i as f64 / n_ticks as f64;
        let yy = oy + ph - frac * ph;
        let tick_val = frac * max_val;
        // gridline
        svg.line(ox, yy, ox + pw, yy, "#e0e0e0", 1.0);
        // tick label
        let label = if tick_val >= 1_000_000.0 {
            format!("{:.1}M", tick_val / 1_000_000.0)
        } else if tick_val >= 1_000.0 {
            format!("{:.1}k", tick_val / 1_000.0)
        } else {
            format!("{:.0}", tick_val)
        };
        svg.text(ox - 4.0, yy + 4.0, &label, 10, "end");
    }

    // Y-axis line
    svg.line(ox, oy, ox, oy + ph, "#666666", 1.0);
    // X-axis line
    svg.line(ox, oy + ph, ox + pw, oy + ph, "#666666", 1.0);

    let n = values.len();
    let bar_w = pw / n as f64 * 0.7;
    let gap = pw / n as f64;

    for (i, (&val, &lbl)) in values.iter().zip(labels.iter()).enumerate() {
        let bh = (val / max_val * ph).max(0.0);
        let bx = ox + i as f64 * gap + (gap - bar_w) / 2.0;
        let by = oy + ph - bh;

        svg.rect(bx, by, bar_w, bh, color, 0.85);

        // Value text inside bar if tall enough
        if bh > 14.0 {
            let val_str = if val >= 1_000_000.0 {
                format!("{:.1}M", val / 1_000_000.0)
            } else if val >= 1_000.0 {
                format!("{:.1}k", val / 1_000.0)
            } else {
                format!("{:.0}", val)
            };
            svg.text(bx + bar_w / 2.0, by + 12.0, &val_str, 9, "middle");
        }

        // X-axis label (truncated to 12 chars)
        let display_lbl = truncate_label(lbl, 12);
        svg.text(bx + bar_w / 2.0, oy + ph + 14.0, &display_lbl, 10, "middle");
    }

    svg.build()
}

/// Scatter plot for UMAP / PCA projections.
pub fn scatter_svg(
    title: &str,
    points: &[(f64, f64)],
    colors: &[&str],
    x_label: &str,
    y_label: &str,
    width: u32,
    height: u32,
) -> String {
    let mut svg = Svg::new(width, height);
    svg.rect(0.0, 0.0, width as f64, height as f64, "#ffffff", 1.0);

    let pw = width as f64 - MARGIN_LEFT - MARGIN_RIGHT;
    let ph = height as f64 - MARGIN_TOP - MARGIN_BOTTOM;
    let ox = MARGIN_LEFT;
    let oy = MARGIN_TOP;

    // Title
    svg.text(ox + pw / 2.0, oy - 10.0, title, 14, "middle");

    // Axis labels
    svg.text(ox + pw / 2.0, oy + ph + 48.0, x_label, 11, "middle");
    // Y-axis label (rotated via transform)
    let ylabel_x = 12.0;
    let ylabel_y = oy + ph / 2.0;
    let escaped_y = escape_svg(y_label);
    svg.elements.push(format!(
        r#"<text x="{ylabel_x:.2}" y="{ylabel_y:.2}" font-size="11" text-anchor="middle" font-family="Arial, sans-serif" transform="rotate(-90,{ylabel_x:.2},{ylabel_y:.2})">{escaped_y}</text>"#
    ));

    if points.is_empty() {
        // Just draw axes
        svg.line(ox, oy, ox, oy + ph, "#666666", 1.0);
        svg.line(ox, oy + ph, ox + pw, oy + ph, "#666666", 1.0);
        return svg.build();
    }

    // Auto-scale with 10% padding
    let (mut xmin, mut xmax, mut ymin, mut ymax) = points.iter().fold(
        (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ),
        |(xn, xx, yn, yx), &(x, y)| (xn.min(x), xx.max(x), yn.min(y), yx.max(y)),
    );
    let xpad = (xmax - xmin).max(1e-9) * 0.1;
    let ypad = (ymax - ymin).max(1e-9) * 0.1;
    xmin -= xpad;
    xmax += xpad;
    ymin -= ypad;
    ymax += ypad;

    let xrange = xmax - xmin;
    let yrange = ymax - ymin;

    // Gridlines
    for i in 0..=4u32 {
        let frac = i as f64 / 4.0;
        let gx = ox + frac * pw;
        let gy = oy + frac * ph;
        svg.line(gx, oy, gx, oy + ph, "#e0e0e0", 1.0);
        svg.line(ox, gy, ox + pw, gy, "#e0e0e0", 1.0);
    }

    // Axes
    svg.line(ox, oy, ox, oy + ph, "#666666", 1.0);
    svg.line(ox, oy + ph, ox + pw, oy + ph, "#666666", 1.0);

    // Axis tick labels
    for i in 0..=4u32 {
        let frac = i as f64 / 4.0;
        let xval = xmin + frac * xrange;
        let yval = ymax - frac * yrange;
        let gx = ox + frac * pw;
        let gy = oy + frac * ph;
        svg.text(gx, oy + ph + 14.0, &format!("{xval:.2}"), 9, "middle");
        svg.text(ox - 4.0, gy + 4.0, &format!("{yval:.2}"), 9, "end");
    }

    // Points
    let fallback_color = "#4e79a7";
    for (i, &(x, y)) in points.iter().enumerate() {
        let px = ox + (x - xmin) / xrange * pw;
        let py = oy + (ymax - y) / yrange * ph;
        let fill = colors.get(i).copied().unwrap_or(fallback_color);
        svg.circle(px, py, 5.0, fill, 0.8);
    }

    svg.build()
}

/// Volcano plot: log₂FC on X, −log₁₀(padj) on Y.
pub fn volcano_svg(title: &str, points: &[(f64, f64, bool)], width: u32, height: u32) -> String {
    let mut svg = Svg::new(width, height);
    svg.rect(0.0, 0.0, width as f64, height as f64, "#ffffff", 1.0);

    let pw = width as f64 - MARGIN_LEFT - MARGIN_RIGHT;
    let ph = height as f64 - MARGIN_TOP - MARGIN_BOTTOM;
    let ox = MARGIN_LEFT;
    let oy = MARGIN_TOP;

    svg.text(ox + pw / 2.0, oy - 10.0, title, 14, "middle");

    // Axis labels
    svg.text(
        ox + pw / 2.0,
        oy + ph + 48.0,
        "log\u{2082} fold change",
        11,
        "middle",
    );
    let ylabel_x = 12.0;
    let ylabel_y = oy + ph / 2.0;
    svg.elements.push(format!(
        "<text x=\"{ylabel_x:.2}\" y=\"{ylabel_y:.2}\" font-size=\"11\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" transform=\"rotate(-90,{ylabel_x:.2},{ylabel_y:.2}\">\u{2212}log\u{2081}\u{2080}(padj)</text>"
    ));

    // Determine axis ranges
    let (xmin, xmax, ymax) = if points.is_empty() {
        (-3.0f64, 3.0f64, 5.0f64)
    } else {
        let xmin = points
            .iter()
            .map(|p| p.0)
            .fold(f64::INFINITY, f64::min)
            .min(-1.6);
        let xmax = points
            .iter()
            .map(|p| p.0)
            .fold(f64::NEG_INFINITY, f64::max)
            .max(1.6);
        let ymax = points
            .iter()
            .map(|p| p.1)
            .fold(f64::NEG_INFINITY, f64::max)
            .max(2.0);
        (xmin, xmax, ymax)
    };
    let ymin = 0.0f64;
    let xrange = xmax - xmin;
    let yrange = ymax - ymin;

    // Gridlines
    for i in 0..=4u32 {
        let frac = i as f64 / 4.0;
        let gy = oy + frac * ph;
        svg.line(ox, gy, ox + pw, gy, "#e0e0e0", 1.0);
    }
    for i in 0..=4u32 {
        let frac = i as f64 / 4.0;
        let gx = ox + frac * pw;
        svg.line(gx, oy, gx, oy + ph, "#e0e0e0", 1.0);
    }

    // Axes
    svg.line(ox, oy, ox, oy + ph, "#666666", 1.0);
    svg.line(ox, oy + ph, ox + pw, oy + ph, "#666666", 1.0);

    // Threshold dashed lines
    let fc_threshold = 1.5_f64;
    let sig_threshold = -0.05_f64.log10(); // ~1.301

    // Vertical FC threshold lines (x = ±1.5)
    for &fc in &[-fc_threshold, fc_threshold] {
        if fc > xmin && fc < xmax {
            let px = ox + (fc - xmin) / xrange * pw;
            // Draw dashed line manually as segments
            let dash_len = 5.0;
            let mut yy = oy;
            let mut draw = true;
            while yy < oy + ph {
                let end = (yy + dash_len).min(oy + ph);
                if draw {
                    svg.line(px, yy, px, end, "#888888", 1.0);
                }
                yy = end;
                draw = !draw;
            }
        }
    }

    // Horizontal significance threshold line (y = 1.301)
    if sig_threshold < ymax {
        let py = oy + (ymax - sig_threshold) / yrange * ph;
        let dash_len = 5.0;
        let mut xx = ox;
        let mut draw = true;
        while xx < ox + pw {
            let end = (xx + dash_len).min(ox + pw);
            if draw {
                svg.line(xx, py, end, py, "#888888", 1.0);
            }
            xx = end;
            draw = !draw;
        }
    }

    // Points
    for &(x, y, sig) in points {
        let px = ox + (x - xmin) / xrange * pw;
        let py = oy + (ymax - y) / yrange * ph;
        let fill = if sig { "#e15759" } else { "#bab0ac" };
        let opacity = if sig { 0.8 } else { 0.4 };
        svg.circle(px, py, 3.0, fill, opacity);
    }

    // Tick labels
    for i in 0..=4u32 {
        let frac = i as f64 / 4.0;
        let xval = xmin + frac * xrange;
        let yval = ymax - frac * yrange;
        svg.text(
            ox + frac * pw,
            oy + ph + 14.0,
            &format!("{xval:.1}"),
            9,
            "middle",
        );
        svg.text(
            ox - 4.0,
            oy + frac * ph + 4.0,
            &format!("{yval:.1}"),
            9,
            "end",
        );
    }

    svg.build()
}

/// Heatmap: blue→white→red color scale, auto-clamped to [-1, 1].
pub fn heatmap_svg(
    title: &str,
    matrix: &[Vec<f64>],
    row_labels: &[&str],
    col_labels: &[&str],
    width: u32,
    height: u32,
) -> String {
    let mut svg = Svg::new(width, height);
    svg.rect(0.0, 0.0, width as f64, height as f64, "#ffffff", 1.0);

    let nrows = matrix.len();
    let ncols = col_labels.len();

    // Title
    svg.text(width as f64 / 2.0, MARGIN_TOP - 10.0, title, 14, "middle");

    if nrows == 0 || ncols == 0 {
        return svg.build();
    }

    // Reserve space for row labels (left) and col labels (top)
    let label_w = 90.0_f64;
    let label_h = 50.0_f64;

    let plot_x = MARGIN_LEFT + label_w;
    let plot_y = MARGIN_TOP + label_h;
    let plot_w = width as f64 - plot_x - MARGIN_RIGHT;
    let plot_h = height as f64 - plot_y - MARGIN_BOTTOM;

    let cell_w = plot_w / ncols as f64;
    let cell_h = plot_h / nrows as f64;

    // Column labels
    for (j, &lbl) in col_labels.iter().enumerate() {
        let cx = plot_x + j as f64 * cell_w + cell_w / 2.0;
        let lbl_short = truncate_label(lbl, 10);
        svg.text(cx, plot_y - 6.0, &lbl_short, 10, "middle");
    }

    // Row labels
    for (i, &lbl) in row_labels.iter().enumerate() {
        let cy = plot_y + i as f64 * cell_h + cell_h / 2.0 + 4.0;
        let lbl_short = truncate_label(lbl, 12);
        svg.text(plot_x - 6.0, cy, &lbl_short, 10, "end");
    }

    // Cells
    for (i, row) in matrix.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            let clamped = val.clamp(-1.0, 1.0);
            let fill = heatmap_color(clamped);
            let cx = plot_x + j as f64 * cell_w;
            let cy = plot_y + i as f64 * cell_h;
            svg.rect(cx, cy, cell_w - 1.0, cell_h - 1.0, &fill, 1.0);

            // Value label
            let text_color = if clamped.abs() > 0.5 {
                "#ffffff"
            } else {
                "#333333"
            };
            svg.elements.push(format!(
                r#"<text x="{:.2}" y="{:.2}" font-size="10" text-anchor="middle" fill="{}" font-family="Arial, sans-serif">{:.2}</text>"#,
                cx + cell_w / 2.0,
                cy + cell_h / 2.0 + 4.0,
                text_color,
                clamped,
            ));
        }
    }

    svg.build()
}

/// Pre-binned histogram for allele frequency distributions.
pub fn histogram_svg(
    title: &str,
    counts: &[u64],
    bin_labels: &[&str],
    color: &str,
    width: u32,
    height: u32,
) -> String {
    let values: Vec<f64> = counts.iter().map(|&c| c as f64).collect();
    bar_chart_svg(title, bin_labels, &values, color, width, height)
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn escape_svg(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn truncate_label(s: &str, max_chars: usize) -> String {
    let chars: Vec<char> = s.chars().collect();
    if chars.len() <= max_chars {
        s.to_string()
    } else {
        let truncated: String = chars[..max_chars - 3].iter().collect();
        format!("{}...", truncated)
    }
}

/// Map a correlation value in [-1, 1] to a CSS hex color.
/// Blue (#4e79a7) at -1, white at 0, red (#e15759) at +1.
fn heatmap_color(r: f64) -> String {
    let r = r.clamp(-1.0, 1.0);
    if r >= 0.0 {
        // white to red
        let t = r;
        let red = 255u8;
        let green = (255.0 * (1.0 - t * 0.12) - t * 88.0).clamp(0.0, 255.0) as u8;
        let blue = (255.0 * (1.0 - t)).clamp(0.0, 255.0) as u8;
        format!("#{:02X}{:02X}{:02X}", red, green, blue)
    } else {
        // blue to white
        let t = -r;
        let blue = 255u8;
        let green = (255.0 * (1.0 - t * 0.14) - t * 40.0).clamp(0.0, 255.0) as u8;
        let red = (255.0 * (1.0 - t * 0.69)).clamp(0.0, 255.0) as u8;
        format!("#{:02X}{:02X}{:02X}", red, green, blue)
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bar_chart_has_svg_tag() {
        let out = bar_chart_svg("T", &["a"], &[1.0], "#4e79a7", 400, 300);
        assert!(out.contains("<svg"), "expected <svg tag");
        assert!(out.contains("</svg>"), "expected </svg> closing tag");
    }

    #[test]
    fn scatter_empty_does_not_panic() {
        let out = scatter_svg("Empty", &[], &["#4e79a7"], "X", "Y", 400, 300);
        assert!(out.contains("<svg"));
    }

    #[test]
    fn heatmap_colors_bounded() {
        // Extreme values should clamp without panic
        let matrix = vec![vec![2.0, -2.0], vec![f64::INFINITY, f64::NEG_INFINITY]];
        let out = heatmap_svg("Test", &matrix, &["r1", "r2"], &["c1", "c2"], 400, 300);
        assert!(out.contains("<svg"));
        // Should not contain any nan/inf in the output
        assert!(!out.contains("nan"), "output should not contain nan");
        assert!(!out.contains("inf"), "output should not contain inf");
    }
}
