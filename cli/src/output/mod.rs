pub mod html;
pub mod json;

pub use html::write_html_report;
pub use json::{build_multiqc_output, write_json, MultiQcOutput};
