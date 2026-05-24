use std::path::Path;

use anyhow::{Context, Result};
use chrono::{DateTime, Utc};

use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use integration_layer::{GeneParadox, GeneRegulatoryProfile, GeneState, Insight, InsightLevel, IntegrationSummary, ParadoxKind};
use proteomics_core::ProteomicsSummary;
use transcriptomics_core::TranscriptomicsSummary;

use crate::context_detect::SampleContext;
use super::circos::generate_circos_svg;
use super::svg::{bar_chart_svg, heatmap_svg, histogram_svg, scatter_svg, volcano_svg};

/// Generate a self-contained HTML report and write it to `{output_dir}/report.html`.
///
/// All charts are rendered as inline SVG — no external dependencies required.
#[allow(clippy::too_many_arguments)]
pub fn write_html_report(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    integration: &IntegrationSummary,
    proteomics: Option<&ProteomicsSummary>,
    sample_context: Option<&SampleContext>,
    generated_at: DateTime<Utc>,
    output_dir: &Path,
) -> Result<()> {
    let html = generate_html(
        genomics,
        transcr,
        epigen,
        integration,
        proteomics,
        sample_context,
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
    sample_context: Option<&SampleContext>,
    generated_at: DateTime<Utc>,
) -> String {
    let head = html_head();
    let summary_cards = html_summary_cards(genomics, transcr, epigen);
    let context_box = sample_context.map(html_sample_context_box).unwrap_or_default();
    let circos_svg = generate_circos_svg(genomics, epigen);
    let variant_chart = html_variant_density_chart(genomics);
    let af_histogram = html_af_histogram(genomics);
    let expression_chart = html_expression_chart(transcr);
    let methylation_chart = html_methylation_chart(epigen);
    let clock_section = html_epigenetic_clock(epigen);
    let volcano_chart = html_volcano_chart(transcr);
    let heatmap = html_correlation_heatmap(integration);
    let pca_chart = html_pca_chart(integration);
    let insights_section = html_insights(&integration.insights);
    let paradoxes_section = html_paradoxes(&integration.paradoxes);
    let cancer_section = html_cancer_genomics(genomics);
    let gene_states_section = html_gene_states(integration);
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
  {context_box}
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
  {clock_section}
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
  {paradoxes_section}
  {cancer_section}
  {gene_states_section}
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
        context_box = context_box,
        circos_svg = circos_svg,
        variant_chart = variant_chart,
        af_histogram = af_histogram,
        expression_chart = expression_chart,
        methylation_chart = methylation_chart,
        clock_section = clock_section,
        volcano_section = if volcano_chart.is_empty() {
            String::new()
        } else {
            format!("<section class=\"section\"><h2>Differential Expression Volcano Plot</h2>{}</section>", volcano_chart)
        },
        heatmap = heatmap,
        pca_chart = pca_chart,
        insights_section = insights_section,
        paradoxes_section = paradoxes_section,
        cancer_section = cancer_section,
        gene_states_section = gene_states_section,
        pathway_table = pathway_table,
        proteomics_section = proteomics_section,
    )
}

/// Render a compact "Sample Context" info box showing auto-detected properties.
fn html_sample_context_box(ctx: &SampleContext) -> String {
    let species = ctx.species.as_str();
    let genomics_assay = ctx
        .genomics_assay
        .as_ref()
        .map(|a| a.as_str())
        .unwrap_or("unknown");
    let epigenomics_assay = ctx
        .epigenomics_assay
        .as_ref()
        .map(|a| a.as_str())
        .unwrap_or("unknown");
    let titv = ctx
        .titv_ratio
        .map(|r| format!("{:.2}", r))
        .unwrap_or_else(|| "n/a".into());
    let burden = ctx
        .somatic_burden_per_mb
        .map(|b| format!("{:.2} /Mb", b))
        .unwrap_or_else(|| "n/a".into());
    let signature = ctx
        .mutation_signature_hint
        .as_deref()
        .unwrap_or("none detected");
    let preset = ctx
        .suggested_preset
        .as_ref()
        .map(|p| format!("{} ({:.0}% confidence)", p.preset, p.confidence * 100.0))
        .unwrap_or_else(|| "none".into());

    let mut warnings_html = String::new();
    for w in ctx.warnings.iter().chain(ctx.concordance.warnings.iter()) {
        warnings_html.push_str(&format!(
            "<div class=\"insight insight-warn\"><span class=\"insight-tag\">[WARN]</span>{}</div>",
            escape_html(w)
        ));
    }

    format!(
        r#"<section class="section">
  <h2>Auto-Detected Sample Context</h2>
  <div class="cards">
    <div class="card">
      <h3>Biology</h3>
      <div class="stat"><span class="stat-label">Species</span><span class="stat-value">{species}</span></div>
      <div class="stat"><span class="stat-label">Genomics assay</span><span class="stat-value">{genomics_assay}</span></div>
      <div class="stat"><span class="stat-label">Epigenomics assay</span><span class="stat-value">{epigenomics_assay}</span></div>
    </div>
    <div class="card">
      <h3>Mutation Statistics</h3>
      <div class="stat"><span class="stat-label">Ti/Tv ratio</span><span class="stat-value">{titv}</span></div>
      <div class="stat"><span class="stat-label">Burden estimate</span><span class="stat-value">{burden}</span></div>
      <div class="stat"><span class="stat-label">Signature hint</span><span class="stat-value">{signature}</span></div>
    </div>
    <div class="card">
      <h3>Suggested Preset</h3>
      <div class="stat"><span class="stat-label">Preset</span><span class="stat-value">{preset}</span></div>
    </div>
  </div>
  {warnings_html}
</section>"#,
        species = species,
        genomics_assay = genomics_assay,
        epigenomics_assay = epigenomics_assay,
        titv = titv,
        burden = burden,
        signature = escape_html(signature),
        preset = escape_html(&preset),
        warnings_html = warnings_html,
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

/// Render the Horvath epigenetic age clock card.
///
/// Returns empty string when no clock result is available.
fn html_epigenetic_clock(e: &EpigenomicsSummary) -> String {
    let ma = match e.methylation_age.as_ref() {
        Some(m) => m,
        None => return String::new(),
    };

    let conf_badge = match ma.confidence.as_str() {
        "HIGH" => format!(
            r#"<span class="badge badge-green">{}</span>"#,
            ma.confidence
        ),
        "MODERATE" => format!(
            r#"<span class="badge badge-yellow">{}</span>"#,
            ma.confidence
        ),
        _ => format!(
            r#"<span class="badge badge-red">{}</span>"#,
            ma.confidence
        ),
    };

    let delta_html = match ma.age_delta {
        Some(d) => {
            let color = if d > 5.0 {
                "var(--red)"
            } else if d < -5.0 {
                "var(--green)"
            } else {
                "var(--text)"
            };
            format!(
                r#"<div class="stat"><span class="stat-label">Age Delta (bio - chron)</span><span class="stat-value" style="color:{color}">{:+.1} years</span></div>"#,
                d
            )
        }
        None => r#"<div class="stat"><span class="stat-label">Age Delta</span><span class="stat-value">N/A (no chronological age)</span></div>"#.to_string(),
    };

    let accel_html = if ma.age_accelerated == Some(true) {
        r#"<div class="insight insight-warn"><span class="insight-tag">[WARN]</span>Epigenetic age acceleration detected — associated with cancer, neurodegeneration, and increased mortality risk.</div>"#
    } else if ma.age_accelerated == Some(false) {
        r#"<div class="insight insight-info"><span class="insight-tag">[INFO]</span>No epigenetic age acceleration detected.</div>"#
    } else {
        ""
    };

    format!(
        r#"<section class="section">
<h2>Epigenetic Age Clock (Horvath 2013)</h2>
<div class="cards">
  <div class="card">
    <h3>Predicted Biological Age</h3>
    <div style="font-size:48px;font-weight:700;color:var(--accent);text-align:center;padding:16px 0">{age:.1}</div>
    <div style="text-align:center;color:#8b949e;font-size:12px">years</div>
  </div>
  <div class="card">
    <h3>Clock Coverage</h3>
    <div class="stat"><span class="stat-label">CpGs found</span><span class="stat-value">{found} / {total}</span></div>
    <div class="stat"><span class="stat-label">Coverage fraction</span><span class="stat-value">{cov:.1}%</span></div>
    <div class="stat"><span class="stat-label">Confidence</span><span class="stat-value">{conf_badge}</span></div>
    {delta_html}
  </div>
</div>
{accel_html}
<p style="font-size:11px;color:#8b949e;margin-top:8px">Based on Horvath 2013 clock (353-CpG elastic net model, 50-site approximation). Coordinates: hg19. Reference: doi:10.1186/gb-2013-14-10-r115.</p>
</section>"#,
        age = ma.biological_age,
        found = ma.cpgs_found,
        total = ma.cpgs_total,
        cov = ma.coverage * 100.0,
        conf_badge = conf_badge,
        delta_html = delta_html,
        accel_html = accel_html,
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

fn html_paradoxes(paradoxes: &[GeneParadox]) -> String {
    if paradoxes.is_empty() {
        return String::new();
    }

    let mut html = String::from(
        r#"<section class="section">
<h2>&#x26A0; Biological Paradoxes</h2>
<table>
<thead><tr>
  <th>Gene</th><th>Type</th><th>Evidence</th><th>Interpretation</th>
</tr></thead>
<tbody>"#,
    );

    for p in paradoxes {
        let (row_style, type_label) = match &p.kind {
            ParadoxKind::MultiHit { n_modalities } => (
                " style=\"background:#2b0f0f\"",
                format!("MultiHit ({n_modalities} modalities)"),
            ),
            ParadoxKind::MethylatedButExpressed => (
                " style=\"background:#2b1a00\"",
                "Methylated + Expressed".to_string(),
            ),
            ParadoxKind::VariantInActiveGene => (
                " style=\"background:#2b2200\"",
                "Variant in Active Gene".to_string(),
            ),
            ParadoxKind::VariantInSilentGene => (
                "",
                "Variant in Silent Gene".to_string(),
            ),
            ParadoxKind::DifferentialWithoutVariant => (
                "",
                "DE without Variant".to_string(),
            ),
            ParadoxKind::VariantWithoutExpression => (
                "",
                "Variant without DE".to_string(),
            ),
        };

        // Build evidence string
        let mut ev_parts: Vec<String> = Vec::new();
        if let Some(tpm) = p.evidence.mean_tpm {
            ev_parts.push(format!("TPM={tpm:.1}"));
        }
        if let Some(meth) = p.evidence.mean_methylation {
            ev_parts.push(format!("Meth={meth:.1}%"));
        }
        if let Some(qual) = p.evidence.max_variant_qual {
            ev_parts.push(format!("QUAL={qual:.0}"));
        }
        if let Some(lfc) = p.evidence.log2_fold_change {
            ev_parts.push(format!("log2FC={lfc:.2}"));
        }
        if let Some(pj) = p.evidence.padj {
            ev_parts.push(format!("padj={pj:.2e}"));
        }
        let evidence_str = ev_parts.join(", ");

        html.push_str(&format!(
            "<tr{row_style}><td><strong>{gene}</strong></td><td>{kind}</td><td>{ev}</td><td>{summary}</td></tr>\n",
            row_style = row_style,
            gene = escape_html(&p.gene),
            kind = escape_html(&type_label),
            ev = escape_html(&evidence_str),
            summary = escape_html(&p.summary),
        ));
    }

    html.push_str("</tbody></table></section>");
    html
}

fn html_gene_states(integration: &IntegrationSummary) -> String {
    let profiles: Vec<&GeneRegulatoryProfile> = integration
        .gene_states
        .iter()
        .filter(|g| g.state != GeneState::Unknown)
        .take(200)
        .collect();

    if profiles.is_empty() {
        return String::new();
    }

    let mut html = String::from(r#"<section class="section">
<h2>&#x1F9EC; Gene Regulatory States</h2>
<p>State assignment based on multi-modal molecular profile (expression &times; methylation &times; genomic variants).
Thresholds: <span style="color:#27ae60">&#9632; Active</span> = TPM &ge; 10 + methylation &lt; 30%;
<span style="color:#7f8c8d">&#9632; Silenced</span> = TPM &lt; 1 + methylation &gt; 70%;
<span style="color:#8e44ad">&#9632; Bivalent</span> = TPM 1&ndash;10 + methylation &gt; 70%;
<span style="color:#f39c12">&#9632; Poised</span> = TPM 1&ndash;10 + methylation &lt; 30%;
<span style="color:#e74c3c">&#9632; Variant-Driven</span> = variant co-occurs with significant DE;
<span style="color:#e67e22">&#9632; Paradoxical</span> = conflicting multi-modal signals.</p>
<table>
<thead><tr>
  <th>Gene</th><th>State</th><th>TPM</th><th>Methylation</th><th>Variant QUAL</th><th>log2FC</th><th>Description</th>
</tr></thead>
<tbody>
"#);

    for g in &profiles {
        let color = g.state.html_color();
        let state_badge = format!(
            r#"<span style="background:{color};color:#fff;padding:2px 7px;border-radius:4px;font-size:0.8em">{}</span>"#,
            g.state.as_str()
        );
        let tpm = g.mean_tpm.map(|v| format!("{v:.1}")).unwrap_or_default();
        let meth = g.mean_methylation.map(|v| format!("{v:.1}%")).unwrap_or_default();
        let qual = g.max_variant_qual.map(|v| format!("{v:.0}")).unwrap_or_default();
        let lfc = g.log2_fold_change.map(|v| format!("{v:+.2}")).unwrap_or_default();

        html.push_str(&format!(
            "<tr><td><strong>{gene}</strong></td><td>{badge}</td><td>{tpm}</td><td>{meth}</td><td>{qual}</td><td>{lfc}</td><td>{desc}</td></tr>\n",
            gene = escape_html(&g.gene),
            badge = state_badge,
            desc = escape_html(&g.description),
        ));
    }

    html.push_str("</tbody></table></section>");
    html
}

fn html_cancer_genomics(g: &GenomicsSummary) -> String {
    // Only render if there is meaningful cancer data
    let has_purity = g.tumor_purity.is_some();
    let has_kataegis = !g.kataegis_loci.is_empty();
    let has_hrd = g.hrd.is_some();
    let has_loh = g.loh_chromosomes.iter().any(|c| c.loh_flagged);

    if !has_purity && !has_kataegis && !has_hrd && !has_loh {
        return String::new();
    }

    let mut html = String::from(
        r#"<section class="section">
<h2>Cancer Genomics Analysis</h2>"#,
    );

    // ── Tumor Purity ──
    if let Some(ref p) = g.tumor_purity {
        let vaf_str = p
            .vaf_purity
            .map(|v| format!("{:.1}%", v * 100.0))
            .unwrap_or_else(|| "N/A".to_string());
        let meth_str = p
            .methylation_purity
            .map(|v| format!("{:.1}%", v * 100.0))
            .unwrap_or_else(|| "N/A".to_string());
        let cons_str = p
            .consensus_purity
            .map(|v| format!("{:.1}%", v * 100.0))
            .unwrap_or_else(|| "N/A".to_string());
        let class_badge = match p.purity_class.as_str() {
            "HIGH" => format!(
                r#"<span class="badge badge-red">{}</span>"#,
                p.purity_class
            ),
            "MODERATE" => format!(
                r#"<span class="badge badge-yellow">{}</span>"#,
                p.purity_class
            ),
            _ => format!(
                r#"<span class="badge badge-green">{}</span>"#,
                p.purity_class
            ),
        };
        let discord_html = if p.discordant {
            r#"<div class="insight insight-warn"><span class="insight-tag">[WARN]</span>Purity estimates discordant between VAF and methylation — possible tumor heterogeneity or technical artifact.</div>"#
        } else {
            ""
        };
        html.push_str(&format!(
            r#"<h3>Tumor Purity Estimation</h3>
<div class="cards">
  <div class="card">
    <h3>Purity Estimates</h3>
    <div class="stat"><span class="stat-label">VAF-based purity</span><span class="stat-value">{vaf}</span></div>
    <div class="stat"><span class="stat-label">Methylation-based purity</span><span class="stat-value">{meth}</span></div>
    <div class="stat"><span class="stat-label">Consensus purity</span><span class="stat-value">{cons}</span></div>
    <div class="stat"><span class="stat-label">Purity class</span><span class="stat-value">{badge}</span></div>
  </div>
</div>
{discord}"#,
            vaf = vaf_str,
            meth = meth_str,
            cons = cons_str,
            badge = class_badge,
            discord = discord_html,
        ));
    }

    // ── Kataegis ──
    html.push_str("<h3>Kataegis Detection</h3>");
    if has_kataegis {
        html.push_str(
            r#"<table><thead><tr>
  <th>Chromosome</th><th>Start</th><th>End</th><th>Mutations</th><th>Mean IMD (bp)</th><th>Dominant Change</th>
</tr></thead><tbody>"#,
        );
        for locus in &g.kataegis_loci {
            html.push_str(&format!(
                "<tr><td>{chrom}</td><td>{start}</td><td>{end}</td><td>{n}</td><td>{imd:.1}</td><td>{change}</td></tr>\n",
                chrom = escape_html(&locus.chrom),
                start = locus.start,
                end = locus.end,
                n = locus.n_mutations,
                imd = locus.geometric_mean_imd,
                change = escape_html(&locus.dominant_change),
            ));
        }
        html.push_str("</tbody></table>");
    } else {
        html.push_str("<p>No kataegis loci detected.</p>");
    }

    // ── HRD ──
    if let Some(ref hrd) = g.hrd {
        let hrd_badge = match hrd.hrd_class.as_str() {
            "HRD-HIGH" => format!(
                r#"<span class="badge badge-red">{}</span>"#,
                hrd.hrd_class
            ),
            "HRD-INTERMEDIATE" => format!(
                r#"<span class="badge badge-yellow">{}</span>"#,
                hrd.hrd_class
            ),
            _ => format!(
                r#"<span class="badge badge-green">{}</span>"#,
                hrd.hrd_class
            ),
        };

        // Inline SVG stacked bar: del_1bp | del_2-5bp | del_6-50bp | ins>3bp
        let bar_w = 400.0_f64;
        let bar_h = 24.0_f64;
        let segments = [
            (hrd.del_1bp_frac, "#4e79a7", "Del 1bp"),
            (hrd.del_2_5bp_frac, "#f28e2b", "Del 2-5bp"),
            (hrd.del_6_50bp_frac, "#e15759", "Del 6-50bp"),
            (hrd.ins_gt3bp_frac, "#76b7b2", "Ins &gt;3bp"),
        ];
        let mut bar_svg =
            format!(r#"<svg viewBox="0 0 {bar_w} {bar_h}" xmlns="http://www.w3.org/2000/svg" style="width:400px;display:block;margin:8px 0">"#, bar_w = bar_w + 120.0, bar_h = bar_h + 4.0);
        let mut x = 0.0_f64;
        for (frac, color, _label) in &segments {
            let w = frac * bar_w;
            if w > 0.0 {
                bar_svg.push_str(&format!(
                    r#"<rect x="{x:.1}" y="0" width="{w:.1}" height="{bar_h}" fill="{color}"/>"#,
                    x = x,
                    w = w,
                    bar_h = bar_h,
                    color = color,
                ));
            }
            x += w;
        }
        bar_svg.push_str("</svg>");
        // Legend
        let legend: String = segments
            .iter()
            .map(|(frac, color, label)| {
                format!(
                    r#"<span style="margin-right:12px"><svg width="12" height="12"><rect width="12" height="12" fill="{color}"/></svg> {label} ({pct:.1}%)</span>"#,
                    color = color,
                    label = label,
                    pct = frac * 100.0,
                )
            })
            .collect();

        let note_html = hrd
            .note
            .as_deref()
            .map(|n| {
                format!(
                    r#"<div class="insight insight-warn"><span class="insight-tag">[NOTE]</span>{}</div>"#,
                    escape_html(n)
                )
            })
            .unwrap_or_default();

        html.push_str(&format!(
            r#"<h3>Homologous Recombination Deficiency (HRD)</h3>
<div class="stat"><span class="stat-label">Total indels</span><span class="stat-value">{indels}</span></div>
<div class="stat"><span class="stat-label">HRD class</span><span class="stat-value">{badge}</span></div>
<div class="stat"><span class="stat-label">HRD-indel score</span><span class="stat-value">{score:.4}</span></div>
{bar_svg}
<div style="font-size:12px;color:#8b949e">{legend}</div>
{note}"#,
            indels = hrd.total_indels,
            badge = hrd_badge,
            score = hrd.hrd_indel_score,
            bar_svg = bar_svg,
            legend = legend,
            note = note_html,
        ));
    }

    // ── LOH ──
    let loh_flagged: Vec<_> = g.loh_chromosomes.iter().filter(|c| c.loh_flagged).collect();
    html.push_str("<h3>Loss of Heterozygosity (LOH)</h3>");
    if loh_flagged.is_empty() {
        html.push_str("<p>No chromosomes with significant LOH detected.</p>");
    } else {
        html.push_str(
            r#"<table><thead><tr>
  <th>Chromosome</th><th>Het Variants</th><th>Median AF Deviation</th><th>Skewed Fraction</th>
</tr></thead><tbody>"#,
        );
        for loh in &loh_flagged {
            html.push_str(&format!(
                "<tr><td>{chrom}</td><td>{n}</td><td>{dev:.4}</td><td>{skew:.3}</td></tr>\n",
                chrom = escape_html(&loh.chrom),
                n = loh.n_het_variants,
                dev = loh.median_af_deviation,
                skew = loh.skewed_fraction,
            ));
        }
        html.push_str("</tbody></table>");
    }

    html.push_str("</section>");
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
