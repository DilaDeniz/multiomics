use std::path::Path;

use anyhow::{Context, Result};
use chrono::{DateTime, Utc};

use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use integration_layer::{Insight, InsightLevel, IntegrationSummary};
use proteomics_core::ProteomicsSummary;
use transcriptomics_core::TranscriptomicsSummary;

use super::circos::generate_circos_svg;
use super::svg::{bar_chart_svg, heatmap_svg, histogram_svg, scatter_svg, volcano_svg};

/// Generate a self-contained HTML report and write it to `{output_dir}/report.html`.
///
/// All charts are rendered as inline SVG — no external dependencies required.
pub fn write_html_report(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    integration: &IntegrationSummary,
    proteomics: Option<&ProteomicsSummary>,
    generated_at: DateTime<Utc>,
    output_dir: &Path,
) -> Result<()> {
    let html = generate_html(
        genomics,
        transcr,
        epigen,
        integration,
        proteomics,
        generated_at,
    );
    let path = output_dir.join("report.html");
    std::fs::write(&path, html)
        .with_context(|| format!("Cannot write HTML report to '{}'", path.display()))?;
    log::info!("HTML report written to '{}'", path.display());
    Ok(())
}

fn generate_html(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    integration: &IntegrationSummary,
    proteomics: Option<&ProteomicsSummary>,
    generated_at: DateTime<Utc>,
) -> String {
    let head = html_head();
    let summary_cards = html_summary_cards(genomics, transcr, epigen);
    let circos_svg = generate_circos_svg(genomics, epigen);
    let variant_chart = html_variant_density_chart(genomics);
    let af_histogram = html_af_histogram(genomics);
    let expression_chart = html_expression_chart(transcr);
    let methylation_chart = html_methylation_chart(epigen);
    let volcano_chart = html_volcano_chart(transcr);
    let heatmap = html_correlation_heatmap(integration);
    let pca_chart = html_pca_chart(integration);
    let insights_section = html_insights(&integration.insights);
    let pathway_table = html_pathway_table(integration);
    let proteomics_section = proteomics.map(html_proteomics_section).unwrap_or_default();

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
{head}
<body>
<div class="container">
  <header>
    <h1>🧬 Multiomics Report</h1>
    <p class="subtitle">Generated: {ts} &nbsp;|&nbsp; Tool version: {ver}</p>
  </header>
  {summary_cards}
  <div class="row-2">
    <section class="section">
      <h2>Genomic Overview</h2>
      {circos_svg}
    </section>
    <section class="section">
      <h2>Variant Density by Chromosome</h2>
      {variant_chart}
    </section>
  </div>
  <section class="section">
    <h2>Allele Frequency Distribution</h2>
    {af_histogram}
  </section>
  <section class="section">
    <h2>Top 20 Expressed Genes</h2>
    {expression_chart}
  </section>
  <section class="section">
    <h2>Per-Chromosome Methylation Profile</h2>
    {methylation_chart}
  </section>
  {volcano_section}
  <div class="row-2">
    <section class="section">
      <h2>Cross-Modality Correlation</h2>
      {heatmap}
    </section>
    <section class="section">
      <h2>PCA Projection (2D)</h2>
      {pca_chart}
    </section>
  </div>
  {insights_section}
  {pathway_table}
  {proteomics_section}
  <footer>
    <p>Multiomics &copy; 2026 — Apache 2.0 License —
       <a href="https://github.com/diladeniz/multiomics">github.com/diladeniz/multiomics</a></p>
  </footer>
</div>
</body>
</html>"#,
        head = head,
        ts = generated_at.format("%Y-%m-%d %H:%M:%S UTC"),
        ver = env!("CARGO_PKG_VERSION"),
        summary_cards = summary_cards,
        circos_svg = circos_svg,
        variant_chart = variant_chart,
        af_histogram = af_histogram,
        expression_chart = expression_chart,
        methylation_chart = methylation_chart,
        volcano_section = if volcano_chart.is_empty() {
            String::new()
        } else {
            format!("<section class=\"section\"><h2>Differential Expression Volcano Plot</h2>{}</section>", volcano_chart)
        },
        heatmap = heatmap,
        pca_chart = pca_chart,
        insights_section = insights_section,
        pathway_table = pathway_table,
        proteomics_section = proteomics_section,
    )
}

fn html_proteomics_section(p: &ProteomicsSummary) -> String {
    let mut rows = String::new();
    for prot in p.top_proteins.iter().take(20) {
        rows.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.1}</td><td>{:.4}</td></tr>\n",
            escape_html(&prot.protein),
            prot.n_psms,
            prot.n_unique_peptides,
            prot.top_score,
            prot.q_value,
        ));
    }

    // Score histogram as a small bar chart.
    let max_bin = p.score_histogram.iter().copied().max().unwrap_or(1).max(1);
    let bars: String = p
        .score_histogram
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            let h = (v as f64 / max_bin as f64 * 80.0) as u32;
            format!(
                "<rect x=\"{x}\" y=\"{y}\" width=\"18\" height=\"{h}\" fill=\"#4e79a7\"/>",
                x = i * 22 + 10,
                y = 90u32.saturating_sub(h),
                h = h,
            )
        })
        .collect();

    format!(
        r#"<section class="section">
  <h2>Proteomics</h2>
  <div class="cards">
    <div class="card"><div class="card-val">{ms2}</div><div class="card-lbl">MS2 Spectra</div></div>
    <div class="card"><div class="card-val">{psms}</div><div class="card-lbl">PSMs (1% FDR)</div></div>
    <div class="card"><div class="card-val">{peps}</div><div class="card-lbl">Peptides (1% FDR)</div></div>
    <div class="card"><div class="card-val">{prots}</div><div class="card-lbl">Proteins (1% FDR)</div></div>
    <div class="card"><div class="card-val">{score:.1}</div><div class="card-lbl">Median Hyperscore</div></div>
  </div>
  <h3>Hyperscore Distribution</h3>
  <svg viewBox="0 0 470 100" xmlns="http://www.w3.org/2000/svg" style="max-width:500px">{bars}</svg>
  <h3>Top Identified Proteins</h3>
  <table>
    <thead><tr><th>Protein</th><th>PSMs</th><th>Peptides</th><th>Top Score</th><th>q-value</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</section>"#,
        ms2 = p.n_ms2,
        psms = p.n_psms_1pct,
        peps = p.n_peptides_1pct,
        prots = p.n_proteins_1pct,
        score = p.median_hyperscore,
        bars = bars,
        rows = rows,
    )
}

fn html_head() -> String {
    r#"<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Multiomics Report</title>
<style>
  :root {
    --bg: #0d1117; --surface: #161b22; --border: #30363d;
    --text: #c9d1d9; --accent: #58a6ff; --green: #3fb950;
    --yellow: #d29922; --red: #f85149; --purple: #bc8cff;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { background: var(--bg); color: var(--text); font-family: 'Segoe UI', sans-serif; font-size: 14px; }
  .container { max-width: 1400px; margin: 0 auto; padding: 24px; }
  header { border-bottom: 1px solid var(--border); padding-bottom: 16px; margin-bottom: 24px; }
  header h1 { font-size: 28px; color: var(--accent); }
  .subtitle { color: #8b949e; margin-top: 4px; }
  .cards { display: flex; gap: 16px; margin-bottom: 24px; flex-wrap: wrap; }
  .card { background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
          padding: 16px; flex: 1; min-width: 280px; }
  .card h3 { color: var(--accent); margin-bottom: 12px; font-size: 14px; text-transform: uppercase; letter-spacing: 0.05em; }
  .stat { display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px solid #21262d; }
  .stat:last-child { border-bottom: none; }
  .stat-label { color: #8b949e; }
  .stat-value { font-weight: 600; color: var(--text); }
  .badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 12px; font-weight: 600; }
  .badge-green { background: #1a4731; color: var(--green); }
  .badge-yellow { background: #3d2b00; color: var(--yellow); }
  .badge-red { background: #3d0f0f; color: var(--red); }
  .section { background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
             padding: 20px; margin-bottom: 16px; }
  .section h2 { font-size: 16px; color: var(--accent); margin-bottom: 16px; }
  .row-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 16px; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
  th { text-align: left; padding: 8px 12px; background: #21262d; color: #8b949e;
       font-weight: 600; text-transform: uppercase; font-size: 11px; letter-spacing: 0.05em; }
  td { padding: 8px 12px; border-bottom: 1px solid #21262d; }
  tr:hover td { background: #21262d; }
  .corr-table td { text-align: center; width: 80px; height: 60px; font-weight: 600; border: 2px solid var(--bg); }
  .insight { padding: 10px 14px; border-radius: 6px; margin-bottom: 8px; border-left: 3px solid; }
  .insight-info { background: #1a2b1a; border-color: var(--green); }
  .insight-warn { background: #2b2200; border-color: var(--yellow); }
  .insight-crit { background: #2b0f0f; border-color: var(--red); }
  .insight-tag { font-weight: 700; margin-right: 8px; font-size: 12px; }
  footer { border-top: 1px solid var(--border); padding-top: 16px; margin-top: 24px;
           color: #8b949e; font-size: 12px; text-align: center; }
  footer a { color: var(--accent); }
  @media (max-width: 800px) { .row-2 { grid-template-columns: 1fr; } }
</style>
</head>"#.to_string()
}

fn html_summary_cards(
    g: &GenomicsSummary,
    t: &TranscriptomicsSummary,
    e: &EpigenomicsSummary,
) -> String {
    let titv_badge = if g.titv_ratio >= 1.8 && g.titv_ratio <= 2.5 {
        format!(
            r#"<span class="badge badge-green">{:.2}</span>"#,
            g.titv_ratio
        )
    } else {
        format!(
            r#"<span class="badge badge-yellow">{:.2}</span>"#,
            g.titv_ratio
        )
    };

    let meth_badge = if e.global_methylation_pct < 40.0 {
        format!(
            r#"<span class="badge badge-red">{:.1}%</span>"#,
            e.global_methylation_pct
        )
    } else {
        format!(
            r#"<span class="badge badge-green">{:.1}%</span>"#,
            e.global_methylation_pct
        )
    };

    let expr_pct = if t.total_genes > 0 {
        t.expressed_genes as f64 / t.total_genes as f64 * 100.0
    } else {
        0.0
    };

    format!(
        r#"<div class="cards">
  <div class="card">
    <h3>🧬 Genomics</h3>
    <div class="stat"><span class="stat-label">Total Variants</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">SNPs</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Indels</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Ti/Tv Ratio</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">High-Impact (QUAL&gt;30)</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Unique Positions</span><span class="stat-value">{}</span></div>
  </div>
  <div class="card">
    <h3>📊 Transcriptomics</h3>
    <div class="stat"><span class="stat-label">Total Genes</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Expressed (TPM≥1)</span><span class="stat-value">{} ({:.1}%)</span></div>
    <div class="stat"><span class="stat-label">Samples</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">DE Genes</span><span class="stat-value">{}</span></div>
  </div>
  <div class="card">
    <h3>🔬 Epigenomics</h3>
    <div class="stat"><span class="stat-label">CpG Sites</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Global Methylation</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">CpG Islands</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Hypermethylated Regions</span><span class="stat-value">{}</span></div>
    <div class="stat"><span class="stat-label">Hypomethylated Regions</span><span class="stat-value">{}</span></div>
  </div>
</div>"#,
        format_num(g.total_variants),
        format_num(g.snp_count),
        format_num(g.indel_count),
        titv_badge,
        format_num(g.high_impact.len() as u64),
        format_num(g.unique_positions),
        format_num(t.total_genes),
        format_num(t.expressed_genes),
        expr_pct,
        t.sample_count,
        t.diff_expr
            .as_ref()
            .map(|d| format_num(d.len() as u64))
            .unwrap_or_else(|| "N/A".to_string()),
        format_num(e.total_sites),
        meth_badge,
        e.cpg_islands.len(),
        e.hypermethylated.len(),
        e.hypomethylated.len(),
    )
}

fn html_variant_density_chart(g: &GenomicsSummary) -> String {
    let mut chroms: Vec<(&String, u64)> = g.per_chrom.iter().map(|(k, v)| (k, v.total)).collect();
    chroms.sort_by(|a, b| chrom_sort_key(a.0).cmp(&chrom_sort_key(b.0)));

    let labels: Vec<&str> = chroms.iter().map(|(c, _)| c.as_str()).collect();
    let values: Vec<f64> = chroms
        .iter()
        .map(|(c, _)| g.per_chrom[*c].total as f64)
        .collect();

    bar_chart_svg(
        "Variant Density by Chromosome",
        &labels,
        &values,
        "#4e79a7",
        900,
        300,
    )
}

fn html_af_histogram(g: &GenomicsSummary) -> String {
    let counts = &g.af_histogram;
    let n = counts.len();
    // Generate bin labels "0.00–0.05", "0.05–0.10", …
    let owned_labels: Vec<String> = (0..n)
        .map(|i| {
            let lo = i as f64 / n as f64;
            let hi = (i + 1) as f64 / n as f64;
            format!("{lo:.2}–{hi:.2}")
        })
        .collect();
    let labels: Vec<&str> = owned_labels.iter().map(|s| s.as_str()).collect();

    histogram_svg(
        "Allele Frequency Distribution",
        counts,
        &labels,
        "#59a14f",
        900,
        300,
    )
}

fn html_expression_chart(t: &TranscriptomicsSummary) -> String {
    let top20: Vec<_> = t.top_100_expressed.iter().take(20).collect();
    let labels: Vec<&str> = top20.iter().map(|(g, _)| g.as_str()).collect();
    let values: Vec<f64> = top20.iter().map(|(_, v)| *v).collect();

    bar_chart_svg(
        "Top 20 Expressed Genes (Mean TPM)",
        &labels,
        &values,
        "#f28e2b",
        900,
        320,
    )
}

fn html_methylation_chart(e: &EpigenomicsSummary) -> String {
    let mut chroms: Vec<_> = e.per_chrom.iter().collect();
    chroms.sort_by(|a, b| chrom_sort_key(a.0).cmp(&chrom_sort_key(b.0)));
    let labels: Vec<&str> = chroms.iter().map(|(c, _)| c.as_str()).collect();
    let values: Vec<f64> = chroms.iter().map(|(_, cm)| cm.mean_methylation).collect();

    bar_chart_svg(
        "Per-Chromosome Methylation Profile",
        &labels,
        &values,
        "#76b7b2",
        900,
        300,
    )
}

/// Emit volcano SVG plot.
///
/// Points are colored: red = significant (padj < 0.05 AND |log2FC| ≥ 1),
/// grey = not significant. Returns empty string when no DE data available.
fn html_volcano_chart(t: &TranscriptomicsSummary) -> String {
    let de = match t.diff_expr.as_ref() {
        Some(d) if !d.is_empty() => d,
        _ => return String::new(),
    };

    // Limit to at most 5000 points for reasonable page size
    let points: Vec<(f64, f64, bool)> = de
        .iter()
        .take(5000)
        .filter(|r| !r.log2_fold_change.is_nan())
        .map(|r| {
            let neg_log10_padj = if r.padj.is_nan() || r.padj <= 0.0 {
                r.log2_fold_change.abs()
            } else {
                -r.padj.log10()
            };
            let sig = !r.padj.is_nan() && r.padj < 0.05 && r.log2_fold_change.abs() >= 1.0;
            (r.log2_fold_change, neg_log10_padj, sig)
        })
        .collect();

    if points.is_empty() {
        return String::new();
    }

    volcano_svg("Differential Expression Volcano Plot", &points, 900, 400)
}

fn html_correlation_heatmap(integration: &IntegrationSummary) -> String {
    let labels = ["Genomics", "Transcriptomics", "Epigenomics"];
    let corr = &integration.correlation_matrix;

    let matrix: Vec<Vec<f64>> = labels
        .iter()
        .enumerate()
        .map(|(i, _)| {
            (0..3)
                .map(|j| corr.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0))
                .collect()
        })
        .collect();

    heatmap_svg(
        "Cross-Modality Correlation",
        &matrix,
        &labels,
        &labels,
        500,
        320,
    )
}

fn html_pca_chart(integration: &IntegrationSummary) -> String {
    let ev = &integration.pca.explained_variance_ratio;
    let ev0 = ev.first().copied().unwrap_or(0.0) * 100.0;
    let ev1 = ev.get(1).copied().unwrap_or(0.0) * 100.0;

    let points: Vec<(f64, f64)> = integration
        .pca
        .points
        .iter()
        .map(|p| (p[0], p[1]))
        .collect();

    let palette = ["#59a14f", "#4e79a7", "#bc8cff"];
    let colors: Vec<&str> = (0..points.len())
        .map(|i| palette.get(i).copied().unwrap_or("#4e79a7"))
        .collect();

    let x_label = format!("PC1 ({ev0:.1}%)");
    let y_label = format!("PC2 ({ev1:.1}%)");

    scatter_svg(
        "PCA Projection (2D)",
        &points,
        &colors,
        &x_label,
        &y_label,
        500,
        320,
    )
}

fn html_insights(insights: &[Insight]) -> String {
    if insights.is_empty() {
        return String::new();
    }
    let mut html = String::from(r#"<section class="section"><h2>💡 Biological Insights</h2>"#);
    for insight in insights {
        let (cls, tag) = match insight.level {
            InsightLevel::Info => ("insight-info", "[INFO]"),
            InsightLevel::Warning => ("insight-warn", "[WARN]"),
            InsightLevel::Critical => ("insight-crit", "[CRIT]"),
        };
        let color = insight.level.color_hex();
        html.push_str(&format!(
            r#"<div class="insight {}"><span class="insight-tag" style="color:{}">{}</span>{}</div>"#,
            cls, color, tag,
            escape_html(&insight.message)
        ));
    }
    html.push_str("</section>");
    html
}

fn html_pathway_table(integration: &IntegrationSummary) -> String {
    if integration.top_pathways.is_empty() {
        return String::new();
    }
    let mut html = String::from(
        r#"<section class="section"><h2>🔗 Pathway Enrichment (Fisher's Exact Test + BH FDR)</h2>
<table><thead><tr>
  <th>Pathway</th><th>Name</th>
  <th>Overlap</th><th>Size</th><th>p-value</th><th>padj</th><th>Score</th>
</tr></thead><tbody>"#,
    );
    for r in integration.top_pathways.iter().take(20) {
        let sig_style = if !r.padj.is_nan() && r.padj < 0.05 {
            " style=\"color:#3fb950;font-weight:600\""
        } else {
            ""
        };
        let pval_str = if r.p_value.is_nan() {
            "N/A".into()
        } else {
            format!("{:.2e}", r.p_value)
        };
        let padj_str = if r.padj.is_nan() {
            "N/A".into()
        } else {
            format!("{:.2e}", r.padj)
        };
        html.push_str(&format!(
            "<tr><td>{}</td><td{}>{}</td><td>{}</td><td>{}</td><td>{}</td><td{}>{}</td><td>{:.4}</td></tr>",
            escape_html(&r.pathway_id),
            sig_style,
            escape_html(&r.pathway_name),
            r.overlap,
            r.pathway_size,
            pval_str,
            sig_style,
            padj_str,
            r.score
        ));
    }
    html.push_str("</tbody></table></section>");
    html
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn format_num(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    let chars: Vec<char> = s.chars().collect();
    for (i, c) in chars.iter().enumerate() {
        if i > 0 && (chars.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(*c);
    }
    result
}

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Sorting key that puts chr1 < chr2 < ... < chr22 < chrX < chrY before lexicographic.
fn chrom_sort_key(c: &str) -> (u32, &str) {
    let stripped = c.strip_prefix("chr").unwrap_or(c);
    match stripped.parse::<u32>() {
        Ok(n) => (n, ""),
        Err(_) => (1000, stripped),
    }
}
