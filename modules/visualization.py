"""
Модуль визуализации: интерактивные карты, термодинамические графики 
и генерация красивых HTML-выравниваний для BLAST.
"""
import plotly.graph_objects as go
import plotly.express as px
from typing import List
from modules.primer_design import PrimerResult

def create_primer_map(
    seq_length: int,
    primers: List[PrimerResult],
    selected_pairs: List[int] = None,
    title: str = "Интерактивная карта праймеров",
) -> go.Figure:
    """Создаёт красивую и детальную интерактивную карту позиций праймеров со стрелками."""
    if selected_pairs is None:
        selected_pairs = list(range(len(primers)))

    fig = go.Figure()

    # Шаблонная последовательность (Матрица отрисована толстой линией)
    fig.add_trace(go.Scatter(
        x=[1, seq_length],
        y=[0, 0],
        mode="lines",
        line=dict(color="#B0BEC5", width=12),
        name="Матрица (Template)",
        hoverinfo="text",
        text=[f"Матрица: 1 - {seq_length} bp", f"Матрица: 1 - {seq_length} bp"],
    ))

    colors = px.colors.qualitative.Pastel + px.colors.qualitative.Set2

    for idx, pair_idx in enumerate(selected_pairs):
        if pair_idx >= len(primers): continue
        p = primers[pair_idx]
        color = colors[idx % len(colors)]
        y_offset = (idx + 1) * 0.8

        fwd_exon = "Да" if getattr(p, "fwd_spans_junction", False) else "Нет"
        rev_exon = "Да" if getattr(p, "rev_spans_junction", False) else "Нет"

        # --- FORWARD ПРАЙМЕР (Стрелка направо 5'->3') ---
        hover_fwd = (
            f"<b>Forward праймер #{pair_idx+1}</b><br>"
            f"Позиция: {p.left_start} ➝ {p.left_start+p.left_length}<br>"
            f"Tm: {p.left_tm:.1f}°C | GC: {p.left_gc:.1f}%<br>"
            f"Стык экзонов: <b>{fwd_exon}</b><br>"
            f"Seq: <span style='font-family: monospace;'>{p.left_seq}</span><br>"
            "<extra></extra>"
        )
        # Тело стрелки
        fig.add_trace(go.Scatter(
            x=[p.left_start, p.left_start + p.left_length], y=[y_offset, y_offset],
            mode="lines", line=dict(color=color, width=6),
            name=f"Fwd {pair_idx+1}", hovertemplate=hover_fwd
        ))
        # Наконечник стрелки
        fig.add_trace(go.Scatter(
            x=[p.left_start + p.left_length], y=[y_offset],
            mode="markers", marker=dict(symbol="triangle-right", size=14, color=color, line=dict(width=1, color="black")),
            showlegend=False, hoverinfo="skip"
        ))

        # --- REVERSE ПРАЙМЕР (Стрелка налево 3'<-5') ---
        hover_rev = (
            f"<b>Reverse праймер #{pair_idx+1}</b><br>"
            f"Позиция: {p.right_start} ➝ {p.right_start-p.right_length}<br>"
            f"Tm: {p.right_tm:.1f}°C | GC: {p.right_gc:.1f}%<br>"
            f"Стык экзонов: <b>{rev_exon}</b><br>"
            f"Seq: <span style='font-family: monospace;'>{p.right_seq}</span><br>"
            "<extra></extra>"
        )
        # Тело стрелки
        fig.add_trace(go.Scatter(
            x=[p.right_start, p.right_start - p.right_length], y=[-y_offset, -y_offset],
            mode="lines", line=dict(color=color, width=6),
            name=f"Rev {pair_idx+1}", hovertemplate=hover_rev
        ))
        # Наконечник стрелки
        fig.add_trace(go.Scatter(
            x=[p.right_start - p.right_length], y=[-y_offset],
            mode="markers", marker=dict(symbol="triangle-left", size=14, color=color, line=dict(width=1, color="black")),
            showlegend=False, hoverinfo="skip"
        ))

        # Область продукта (полупрозрачный фон)
        fig.add_vrect(
            x0=p.left_start,
            x1=p.right_start,
            fillcolor=color,
            opacity=0.15,
            line_width=1.5,
            line_dash="dot",
            annotation_text=f"<b>Продукт {p.product_size} bp</b>",
            annotation_position="top",
            annotation_font_size=12,
            annotation_font_color=color
        )

    fig.update_layout(
        title=dict(text=title, font=dict(size=22, color="#1A73E8")),
        xaxis=dict(title="Позиция в гене (п.н.)", showgrid=True, gridcolor="#ECEFF1", zeroline=False),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=True, zerolinecolor="#90A4AE", zerolinewidth=3),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.05, xanchor="right", x=1),
        height=400 + len(selected_pairs) * 50,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=20, r=20, t=100, b=20),
        hovermode="closest",
        clickmode="event+select"
    )

    return fig


def add_exons_to_map(fig, exons_1b):
    """Красиво отрисовывает экзоны на карте матрицы с чередующимися цветами."""
    bg_colors = ["#e0f2f1", "#b2dfdb"] 
    for idx, (s1, e1) in enumerate(exons_1b):
        c = bg_colors[idx % 2]
        fig.add_vrect(
            x0=s1, x1=e1, 
            fillcolor=c, opacity=0.7, 
            line_width=1, line_color="#00695c", 
            annotation_text=f"<b>Экзон {idx+1}</b>", annotation_position="top left",
            annotation_font_size=11, annotation_font_color="#004d40", layer="below"
        )
    
    # Красная пунктирная линия, точно указывающая на стык (Junction)
    junctions = [e for (s, e) in exons_1b[:-1]] if exons_1b else []
    for j in junctions:
        fig.add_vline(
            x=j, line_width=2.5, line_dash="dash", line_color="#d32f2f",
            annotation_text="Стык", annotation_position="bottom right",
            annotation_font_color="#d32f2f", layer="above"
        )
    return fig


def create_tm_comparison_chart(primers: List[PrimerResult]) -> go.Figure:
    pairs = [f"#{p.pair_id+1}" for p in primers]
    fig = go.Figure()
    
    fig.add_trace(go.Bar(name="Forward Tm", x=pairs, y=[p.left_tm for p in primers], marker_color="#4285F4", text=[f"{p.left_tm:.1f}°" for p in primers], textposition="auto", marker_line_color="#1967D2", marker_line_width=1))
    fig.add_trace(go.Bar(name="Reverse Tm", x=pairs, y=[p.right_tm for p in primers], marker_color="#EA4335", text=[f"{p.right_tm:.1f}°" for p in primers], textposition="auto", marker_line_color="#C5221F", marker_line_width=1))

    fig.update_layout(
        title="🌡️ Температура плавления (Tm)",
        yaxis_title="Tm (°C)",
        barmode="group",
        plot_bgcolor="white",
        height=350,
        yaxis=dict(range=[min([p.left_tm for p in primers]+[p.right_tm for p in primers])-5, max([p.left_tm for p in primers]+[p.right_tm for p in primers])+5])
    )
    return fig


def create_gc_comparison_chart(primers: List[PrimerResult]) -> go.Figure:
    pairs = [f"#{p.pair_id+1}" for p in primers]
    fig = go.Figure()
    
    fig.add_trace(go.Bar(name="Forward GC%", x=pairs, y=[p.left_gc for p in primers], marker_color="#34A853", text=[f"{p.left_gc:.1f}%" for p in primers], textposition="auto"))
    fig.add_trace(go.Bar(name="Reverse GC%", x=pairs, y=[p.right_gc for p in primers], marker_color="#FBBC05", text=[f"{p.right_gc:.1f}%" for p in primers], textposition="auto"))

    fig.add_hrect(y0=45, y1=60, fillcolor="#34A853", opacity=0.1, line_width=0, annotation_text="RT-qPCR Optimum (45-60%)", annotation_position="top left")

    fig.update_layout(
        title="🧪 GC-содержание",
        yaxis_title="GC %",
        barmode="group",
        plot_bgcolor="white",
        height=350,
        yaxis=dict(range=[20, 80])
    )
    return fig


def create_penalty_chart(primers: List[PrimerResult]) -> go.Figure:
    pairs = [f"#{p.pair_id+1}" for p in primers]
    penalties = [p.penalty for p in primers]
    colors = ["#34A853" if pen < 1 else "#FBBC05" if pen < 3 else "#EA4335" for pen in penalties]

    fig = go.Figure(go.Bar(x=pairs, y=penalties, marker_color=colors, text=[f"{pen:.2f}" for pen in penalties], textposition="auto"))
    fig.update_layout(
        title="⚖️ Штраф Primer3 (Меньше = Лучше)",
        yaxis_title="Penalty",
        plot_bgcolor="white",
        height=300,
    )
    return fig


def render_blast_alignment_html(query_seq: str, subject_seq: str, q_start: int, s_start: int, match_seq: str = None) -> str:
    """
    Генерирует красивый HTML блок для выравнивания BLAST.
    Мисматчи ярко подсвечиваются красным. 
    """
    if not query_seq or not subject_seq:
        return "<span style='color:#F44336;'>Нет данных о выравнивании.</span>"
        
    html = '<div style="background-color: #1E1E1E; color: #D4D4D4; padding: 15px; border-radius: 8px; font-family: \'Courier New\', Courier, monospace; overflow-x: auto; font-size: 14px; line-height: 1.5; box-shadow: inset 0 0 10px rgba(0,0,0,0.5);">'
    
    q_line = f"<span style='color:#569CD6;'>Query  {q_start:<4}</span> "
    m_line = f"            "
    s_line = f"<span style='color:#4EC9B0;'>Sbjct  {s_start:<4}</span> "

    mismatches_count = 0

    for i, (q, s) in enumerate(zip(query_seq, subject_seq)):
        m_char = match_seq[i] if match_seq and i < len(match_seq) else ('|' if q == s else ' ')
        if q == s:
            q_line += f"<span style='color:#D4D4D4;'>{q}</span>"
            m_line += f"<span style='color:#608B4E;'>{m_char}</span>"
            s_line += f"<span style='color:#D4D4D4;'>{s}</span>"
        else:
            mismatches_count += 1
            # Красная подсветка для мисматчей/гэпов
            q_char = q if q != ' ' else '-'
            s_char = s if s != ' ' else '-'
            q_line += f"<span style='background-color:#E51400; color:#FFFFFF; font-weight:bold; padding:0 1px; border-radius:2px;'>{q_char}</span>"
            m_line += " "
            s_line += f"<span style='background-color:#E51400; color:#FFFFFF; font-weight:bold; padding:0 1px; border-radius:2px;'>{s_char}</span>"

    html += f"{q_line}<br>{m_line}<br>{s_line}"
    html += "</div>"
    
    # Добавляем плашку со статусом
    if mismatches_count == 0:
        badge = "<div style='margin-bottom:8px;'><span style='background-color:#4CAF50; color:white; padding:4px 8px; border-radius:4px; font-size:12px; font-weight:bold;'>100% Идентичность</span></div>"
    else:
        badge = f"<div style='margin-bottom:8px;'><span style='background-color:#F44336; color:white; padding:4px 8px; border-radius:4px; font-size:12px; font-weight:bold;'>Мисматчей: {mismatches_count}</span></div>"

    return badge + html