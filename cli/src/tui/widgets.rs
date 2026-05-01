use ratatui::{
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Gauge, List, ListItem, Paragraph, Wrap},
    Frame,
};

use super::app::{AppState, Phase};

/// Render the full BioMultiOmics TUI to the given frame.
pub fn render(frame: &mut Frame, state: &AppState) {
    let area = frame.size();

    // Outer block
    let outer = Block::default()
        .title(Span::styled(
            " BioMultiOmics ",
            Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD),
        ))
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Cyan));
    let inner = outer.inner(area);
    frame.render_widget(outer, area);

    // Split into: top status bar, middle content, bottom key hints
    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1), // status bar
            Constraint::Min(10),   // main content
            Constraint::Length(1), // key hints
        ])
        .split(inner);

    render_status_bar(frame, rows[0], state);
    render_main_content(frame, rows[1], state);
    render_key_hints(frame, rows[2]);
}

fn render_status_bar(frame: &mut Frame, area: Rect, state: &AppState) {
    let elapsed = format_duration(state.elapsed_secs);
    let eta = state
        .eta_secs
        .map(|s| format!("ETA: {}", format_duration(s)))
        .unwrap_or_else(|| "ETA: --:--".to_string());

    let text = Line::from(vec![
        Span::styled(
            format!("  Phase: {:30}", state.phase.label()),
            Style::default().fg(Color::White),
        ),
        Span::styled(
            format!("Elapsed: {}  {}", elapsed, eta),
            Style::default().fg(Color::DarkGray),
        ),
    ]);
    frame.render_widget(Paragraph::new(text), area);
}

fn render_main_content(frame: &mut Frame, area: Rect, state: &AppState) {
    // Split horizontally: left gauges (2/3), right insights (1/3)
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(65), Constraint::Percentage(35)])
        .split(area);

    render_progress_gauges(frame, cols[0], state);
    render_insights_panel(frame, cols[1], state);
}

fn render_progress_gauges(frame: &mut Frame, area: Rect, state: &AppState) {
    let block = Block::default()
        .title(" Progress ")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::DarkGray));
    let inner = block.inner(area);
    frame.render_widget(block, area);

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Min(0),
        ])
        .split(inner);

    render_gauge(
        frame,
        rows[0],
        "GENOMICS",
        state.genomics_pct,
        state.genomics_rps,
        Color::Green,
        matches!(state.phase, Phase::Genomics),
    );
    render_gauge(
        frame,
        rows[1],
        "TRANSCRIPTOMICS",
        state.transcr_pct,
        state.transcr_rps,
        Color::Blue,
        matches!(state.phase, Phase::Transcriptomics),
    );
    render_gauge(
        frame,
        rows[2],
        "EPIGENOMICS",
        state.epigen_pct,
        state.epigen_rps,
        Color::Magenta,
        matches!(state.phase, Phase::Epigenomics),
    );
    render_gauge(
        frame,
        rows[3],
        "INTEGRATION",
        state.integration_pct,
        0.0,
        Color::Yellow,
        matches!(state.phase, Phase::Integration),
    );
}

fn render_gauge(
    frame: &mut Frame,
    area: Rect,
    label: &str,
    pct: f64,
    rps: f64,
    color: Color,
    active: bool,
) {
    let ratio = (pct / 100.0).clamp(0.0, 1.0);
    let rps_str = if rps > 0.0 {
        format!("{:.0} rec/s", rps)
    } else {
        "waiting...".to_string()
    };
    let gauge_label = format!("{:<18} {:>12}", label, rps_str);

    let style = if active {
        Style::default().fg(color).add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(color)
    };

    let gauge = Gauge::default()
        .block(Block::default())
        .gauge_style(style)
        .ratio(ratio)
        .label(gauge_label);
    frame.render_widget(gauge, area);
}

fn render_insights_panel(frame: &mut Frame, area: Rect, state: &AppState) {
    let block = Block::default()
        .title(" Live Insights ")
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::DarkGray));
    let inner = block.inner(area);
    frame.render_widget(block, area);

    let items: Vec<ListItem> = state
        .insights_live
        .iter()
        .map(|msg| {
            let color = if msg.starts_with("[CRIT]") {
                Color::Red
            } else if msg.starts_with("[WARN]") {
                Color::Yellow
            } else {
                Color::Green
            };
            let line = Line::from(Span::styled(msg.clone(), Style::default().fg(color)));
            ListItem::new(line)
        })
        .collect();

    let list = List::new(items);
    frame.render_widget(
        Paragraph::new(
            state
                .insights_live
                .iter()
                .map(|msg| {
                    let color = if msg.starts_with("[CRIT]") {
                        Color::Red
                    } else if msg.starts_with("[WARN]") {
                        Color::Yellow
                    } else {
                        Color::Green
                    };
                    Line::from(Span::styled(msg.clone(), Style::default().fg(color)))
                })
                .collect::<Vec<_>>(),
        )
        .wrap(Wrap { trim: true }),
        inner,
    );

    // Suppress unused variable warning for list
    drop(list);
}

fn render_key_hints(frame: &mut Frame, area: Rect) {
    let hints = Paragraph::new(Line::from(vec![
        Span::styled("  q", Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
        Span::styled(": quit  ", Style::default().fg(Color::DarkGray)),
        Span::styled("p", Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
        Span::styled(": pause  ", Style::default().fg(Color::DarkGray)),
        Span::styled("Ctrl-C", Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
        Span::styled(": exit", Style::default().fg(Color::DarkGray)),
    ]))
    .alignment(Alignment::Left);
    frame.render_widget(hints, area);
}

fn format_duration(secs: u64) -> String {
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    if h > 0 {
        format!("{:02}:{:02}:{:02}", h, m, s)
    } else {
        format!("{:02}:{:02}", m, s)
    }
}
