"""HTML reporting for KmerSutra."""

from __future__ import annotations

import html
from pathlib import Path


def _table_html(*, title: str, records: list[dict[str, object]]) -> str:
    """Build a simple HTML table.

    Parameters
    ----------
    title : str
        Table title.
    records : list[dict[str, object]]
        Records to render.

    Returns
    -------
    str
        HTML fragment.
    """
    if not records:
        return f"<h2>{html.escape(title)}</h2><p>No records available.</p>"
    columns = list(records[0].keys())
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    rows = []
    for record in records:
        cells = "".join(
            f"<td>{html.escape(str(record.get(column, '')))}</td>"
            for column in columns
        )
        rows.append(f"<tr>{cells}</tr>")
    return (
        f"<h2>{html.escape(title)}</h2>"
        "<div class='table-wrap'><table>"
        f"<thead><tr>{header}</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def build_html_report(
    *,
    title: str,
    panel_summary: list[dict[str, object]] | None = None,
    hit_summary: list[dict[str, object]] | None = None,
    detection_calls: list[dict[str, object]] | None = None,
) -> str:
    """Build a standalone KmerSutra HTML report.

    Parameters
    ----------
    title : str
        Report title.
    panel_summary : list[dict[str, object]] | None, optional
        Panel summary records.
    hit_summary : list[dict[str, object]] | None, optional
        Hit summary records.
    detection_calls : list[dict[str, object]] | None, optional
        Detection-call records.

    Returns
    -------
    str
        Complete HTML report.
    """
    style = """
    <style>
      body { font-family: Arial, Helvetica, sans-serif; margin: 0; background: #f7f9fc; color: #1a1a1a; }
      .container { max-width: 1400px; margin: 0 auto; padding: 32px; background: #ffffff; }
      h1, h2 { color: #1f4e79; }
      .note { background: #eef5fb; border-left: 5px solid #1f4e79; padding: 12px 16px; margin: 18px 0; }
      .table-wrap { overflow-x: auto; border: 1px solid #dbe5f0; border-radius: 8px; margin-bottom: 28px; }
      table { border-collapse: collapse; width: 100%; font-size: 13px; }
      th { background: #1f4e79; color: white; padding: 8px; text-align: left; position: sticky; top: 0; }
      td { border-top: 1px solid #e4edf5; padding: 7px; vertical-align: top; }
      tr:nth-child(even) { background: #f9fbfd; }
    </style>
    """
    sections = [
        _table_html(title="Panel summary", records=panel_summary or []),
        _table_html(title="Hit summary", records=hit_summary or []),
        _table_html(title="Detection calls", records=detection_calls or []),
    ]
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title>{style}</head><body>"
        "<div class='container'>"
        f"<h1>{html.escape(title)}</h1>"
        "<p><strong>KmerSutra</strong> diagnostic k-mer report.</p>"
        "<div class='note'>Confidence scores are heuristic until calibrated "
        "against spike-in truth data. High-confidence calls require sufficient "
        "species-specific evidence and low conflicting signal.</div>"
        + "".join(sections)
        + "</div></body></html>"
    )


def write_html_report(
    *,
    output_path: str | Path,
    title: str,
    panel_summary: list[dict[str, object]] | None = None,
    hit_summary: list[dict[str, object]] | None = None,
    detection_calls: list[dict[str, object]] | None = None,
) -> None:
    """Write a KmerSutra HTML report.

    Parameters
    ----------
    output_path : str or pathlib.Path
        Output HTML path.
    title : str
        Report title.
    panel_summary : list[dict[str, object]] | None, optional
        Panel summary records.
    hit_summary : list[dict[str, object]] | None, optional
        Hit summary records.
    detection_calls : list[dict[str, object]] | None, optional
        Detection-call records.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        build_html_report(
            title=title,
            panel_summary=panel_summary,
            hit_summary=hit_summary,
            detection_calls=detection_calls,
        ),
        encoding="utf-8",
    )
