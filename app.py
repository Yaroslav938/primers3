"""
PrimerDesigner — Streamlit-приложение для подбора и проверки праймеров.
Полная версия + Комплементарность в таблице!
"""

import streamlit as st
import pandas as pd
import sys
import os
import re

# Добавляем корневую директорию проекта в path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from modules.primer_design import (
    design_primers, DesignParams, PrimerResult,
    get_qpcr_defaults, get_standard_defaults, get_design_errors,
)
from modules.blast_check import (
    check_primer_pair_specificity, ORGANISM_GROUPS, BLAST_DATABASES, format_blast_summary,
)
from modules.utils import (
    validate_sequence, parse_fasta, primers_to_dataframe, export_to_csv, generate_primer_report,
    calc_tm, calc_gc, calc_hairpin_tm, calc_homodimer_tm, calc_heterodimer_tm, reverse_complement,
    run_virtual_pcr
)
from modules.visualization import (
    create_primer_map, create_tm_comparison_chart, create_gc_comparison_chart, create_penalty_chart,
)

# Пытаемся импортировать NCBI Fetcher
try:
    from modules.ncbi_fetcher import fetch_sequence_and_exons
    HAS_NCBI = True
except ImportError:
    HAS_NCBI = False


# ─── Настройка страницы ────────────────────────────────────────
st.set_page_config(
    page_title="🧬 PrimerDesigner Pro", 
    page_icon="🧬", 
    layout="wide", 
    initial_sidebar_state="expanded"
)

# ─── Вспомогательные функции ──────────────────────────────────
def parse_exons_1based(text: str):
    exons = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line: continue
        m = re.match(r"^(\d+)\s*[-:.]{1,2}\s*(\d+)$", line)
        if not m: raise ValueError(f"Неверный формат экзона: '{line}'. Ожидаю 'start-end' (например, 1-120).")
        s, e = int(m.group(1)), int(m.group(2))
        if s < 1 or e < 1 or e < s: raise ValueError(f"Неверный диапазон экзона: '{line}'.")
        exons.append((s, e))
    exons.sort(key=lambda x: (x[0], x[1]))
    return exons

def exon_junctions_0based(exons_1b):
    if not exons_1b: return []
    return [e for (s, e) in exons_1b[:-1]]

def spans_junction(primer_start: int, primer_len: int, junction: int, min_tail: int) -> bool:
    a, b = primer_start, primer_start + primer_len
    left_tail, right_tail = junction - a, b - junction
    return (a < junction < b) and (left_tail >= min_tail) and (right_tail >= min_tail)

def annotate_primers_with_junctions(primers, junctions, min_tail: int):
    out = []
    for p in primers:
        d = p.__dict__.copy()
        rev_start, rev_len = d["right_start"] - d["right_length"], d["right_length"]
        d["fwd_spans_junction"] = bool([j for j in junctions if spans_junction(d["left_start"], d["left_length"], j, min_tail)])
        d["rev_spans_junction"] = bool([j for j in junctions if spans_junction(rev_start, rev_len, j, min_tail)])
        out.append(d)
    return out

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
    for j in exon_junctions_0based(exons_1b):
        fig.add_vline(
            x=j, line_width=2.5, line_dash="dash", line_color="#d32f2f",
            annotation_text="Стык", annotation_position="bottom right",
            annotation_font_color="#d32f2f", layer="above"
        )
    return fig

def render_blast_alignment_html(query_seq: str, subject_seq: str, q_start: int, s_start: int, match_seq: str = None) -> str:
    """Генерирует HTML-блок графического выравнивания с подсветкой ошибок."""
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
            q_char = q if q != ' ' else '-'
            s_char = s if s != ' ' else '-'
            q_line += f"<span style='background-color:#E51400; color:#FFFFFF; font-weight:bold; padding:0 1px; border-radius:2px;'>{q_char}</span>"
            m_line += " "
            s_line += f"<span style='background-color:#E51400; color:#FFFFFF; font-weight:bold; padding:0 1px; border-radius:2px;'>{s_char}</span>"

    html += f"{q_line}<br>{m_line}<br>{s_line}</div>"
    
    if mismatches_count == 0:
        badge = "<div style='margin-bottom:8px;'><span style='background-color:#4CAF50; color:white; padding:4px 8px; border-radius:4px; font-size:12px; font-weight:bold;'>100% Идентичность</span></div>"
    else:
        badge = f"<div style='margin-bottom:8px;'><span style='background-color:#F44336; color:white; padding:4px 8px; border-radius:4px; font-size:12px; font-weight:bold;'>Мисматчей: {mismatches_count}</span></div>"

    return badge + html

def parse_region(text: str):
    """Парсит текст региона вида 'start, length'."""
    if text:
        m = re.match(r"(\d+)[,\s]+(\d+)", text.strip())
        if m: 
            return int(m.group(1)), int(m.group(2))
    return None, None


# ─── CSS-стилизация ────────────────────────────────────────────
st.markdown("""
<style>
    .main-header { font-size: 2.5rem; color: #1a73e8; font-weight: 800; letter-spacing: -0.5px; margin-bottom: 0px; }
    .sub-header { font-size: 1.1rem; color: #5F6368; margin-bottom: 25px; font-weight: 500; }
    .primer-box { background-color: #f8f9fa; border-left: 4px solid #1a73e8; padding: 10px 15px; margin-bottom: 10px; border-radius: 4px; font-family: monospace; font-size: 1.1rem;}
    .stat-box { background: #E8F0FE; padding: 15px; border-radius: 10px; border-left: 5px solid #1a73e8; margin-bottom: 10px; }
    div[data-testid="stMetricValue"] { font-size: 1.8rem; font-weight: 700; color: #202124; }
    div[data-testid="stMetricLabel"] { font-size: 0.95rem; color: #5F6368; font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px;}
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">🧬 PrimerDesigner Pro</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Профессиональная среда подбора праймеров, BLAST-выравнивания и In-Silico PCR</div>', unsafe_allow_html=True)


# ─── Инициализация session_state ──────────────────────────────
if "primers" not in st.session_state: st.session_state.primers = []
if "primers_df" not in st.session_state: st.session_state.primers_df = None
if "blast_results" not in st.session_state: st.session_state.blast_results = {}
if "sequence" not in st.session_state: st.session_state.sequence = ""
if "exons_1b" not in st.session_state: st.session_state.exons_1b = []
if "loaded_seq" not in st.session_state: st.session_state.loaded_seq = ""
if "seq_name" not in st.session_state: st.session_state.seq_name = "input_sequence"


# ─── Боковая панель ───────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Параметры дизайна")
    mode = st.radio("Режим:", ["Стандартный ПЦР", "qPCR по инструкции (RT-qPCR)"], index=1)
    defaults = get_qpcr_defaults() if mode == "qPCR по инструкции (RT-qPCR)" else get_standard_defaults()

    with st.expander("📐 Длина праймера и продукта", expanded=True):
        col1, col2, col3 = st.columns(3)
        with col1: primer_min_size = st.number_input("Мин", 12, 35, int(defaults.primer_min_size))
        with col2: primer_opt_size = st.number_input("Опт", 12, 35, int(defaults.primer_opt_size))
        with col3: primer_max_size = st.number_input("Макс", 12, 35, int(defaults.primer_max_size))
        
        c1, c2 = st.columns(2)
        with c1: product_min = st.number_input("Прод. Мин", 40, 5000, int(defaults.product_size_min))
        with c2: product_max = st.number_input("Прод. Макс", 40, 5000, int(defaults.product_size_max))

    with st.expander("🌡️ Термодинамика (Tm & GC)", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1: primer_min_tm = st.number_input("Мин Tm", 40.0, 80.0, float(defaults.primer_min_tm), 0.5)
        with col2: primer_opt_tm = st.number_input("Опт Tm", 40.0, 80.0, float(defaults.primer_opt_tm), 0.5)
        with col3: primer_max_tm = st.number_input("Макс Tm", 40.0, 80.0, float(defaults.primer_max_tm), 0.5)
        primer_max_tm_diff = st.number_input("Макс. разница Tm пары (ΔTm)", 0.0, 10.0, float(getattr(defaults, 'primer_max_tm_diff', 2.0)), 0.5)

        gc_range = st.slider("Диапазон GC%", 20.0, 80.0, (float(defaults.primer_min_gc), float(defaults.primer_max_gc)), 1.0)
        primer_gc_clamp = st.number_input("GC-clamp (нуклеотидов на 3'-конце)", 0, 5, int(getattr(defaults, 'primer_gc_clamp', 0)))

    with st.expander("🧪 Концентрации и Соли (Влияют на Tm)", expanded=False):
        st.caption("Позволяет рассчитать точную Tm для вашей ПЦР-смеси.")
        primer_dna_conc = st.number_input("Конц. праймеров (nM)", 10.0, 1000.0, float(defaults.primer_dna_conc), step=10.0)
        primer_salt_mono = st.number_input("Моновалентные катионы K+/Na+ (mM)", 0.0, 200.0, float(defaults.primer_salt_monovalent))
        primer_salt_div = st.number_input("Дивалентные катионы Mg2+ (mM)", 0.0, 15.0, float(defaults.primer_salt_divalent), step=0.1)
        primer_dntp_conc = st.number_input("Конц. dNTPs (mM)", 0.0, 5.0, float(defaults.primer_dntp_conc), step=0.1)

    with st.expander("⚙️ Вторичные структуры и повторы", expanded=False):
        st.caption("Ограничение образования шпилек и димеров.")
        max_self_any = st.number_input("Max Self Any (Общая компл.)", 0.0, 20.0, float(defaults.primer_max_self_any))
        max_self_end = st.number_input("Max 3' Self (3'-Димеры)", 0.0, 20.0, float(defaults.primer_max_self_end))
        max_poly_x = st.number_input("Poly-X (Макс. повторов одного нуклеотида)", 1, 10, int(defaults.primer_max_poly_x))

    with st.expander("🧬 Экзоны (Exon-exon junction)", expanded=(mode == "qPCR по инструкции (RT-qPCR)")):
        enforce_junction = st.checkbox("Требовать праймер на стыке экзонов", value=False, disabled=(mode != "qPCR по инструкции (RT-qPCR)"))
        min_tail = st.number_input("Мин. хвост на экзоне (п.н.)", 1, 12, 4, disabled=(not enforce_junction))

    num_return = st.number_input("Количество возвращаемых пар", 1, 50, int(defaults.num_return))

    st.divider()
    if st.button("🗑️ Сбросить результаты", use_container_width=True):
        for k in ["primers", "primers_df", "blast_results", "sequence", "exons_1b", "loaded_seq", "seq_name"]:
            if k in st.session_state: del st.session_state[k]
        st.rerun()


# ─── Вкладки ────────────────────────────────────────────────────
tab_design, tab_blast, tab_manual, tab_results = st.tabs([
    "📐 Подбор", "🎯 BLAST Выравнивание", "🛠️ In-Silico PCR & Проверка", "📊 Экспорт и Графики"
])


# ─── Вкладка 1: Подбор праймеров ───────────────────────────────
with tab_design:
    st.header("Ввод последовательности")
    
    input_options = ["Вставить текст", "Загрузить FASTA"]
    if HAS_NCBI:
        input_options.insert(0, "Скачать по NCBI Accession (АВТО)")
        
    input_method = st.radio("Способ получения:", input_options, horizontal=True)

    if input_method == "Скачать по NCBI Accession (АВТО)":
        st.info("💡 Автоматически скачает последовательность мРНК и координаты экзонов из GenBank!")
        col_id, col_btn = st.columns([3, 1])
        with col_id:
            accession_id = st.text_input("Accession ID (например, NM_212832.2)", placeholder="NM_212832.2", label_visibility="collapsed")
        with col_btn:
            if st.button("📥 Загрузить", type="primary", use_container_width=True):
                with st.spinner("Связываемся с NCBI..."):
                    seq, exons, err = fetch_sequence_and_exons(accession_id)
                    if err:
                        st.error(err)
                    else:
                        st.session_state.loaded_seq = seq
                        st.session_state.exons_1b = exons
                        st.session_state.seq_name = accession_id
                        st.success(f"✅ Успешно! Длина: {len(seq)} п.н. Найдено экзонов: {len(exons)}")
        
        seq_input = st.session_state.get('loaded_seq', '')
        
    elif input_method == "Загрузить FASTA":
        uploaded_file = st.file_uploader("Файл (.fasta, .txt)", type=["fasta", "fas", "fa", "txt"])
        if uploaded_file is not None:
            st.session_state.loaded_seq = uploaded_file.getvalue().decode("utf-8")
            st.session_state.seq_name = uploaded_file.name
        seq_input = st.session_state.get('loaded_seq', '')
        
    else: # Вставить текст
        seq_input = st.text_area("Вставьте последовательность", value=st.session_state.get('loaded_seq', ''), height=200)

    # Блок работы с регионами и экзонами
    with st.expander("📍 Регионы поиска и Экзоны (опционально)", expanded=False):
        st.markdown("**Настройка регионов (Target / Excluded Region)**")
        st.caption("Формат: `старт, длина` (например `100, 50`).")
        rc1, rc2 = st.columns(2)
        target_reg = rc1.text_input("Включить регион (Target)", placeholder="Например: 50, 100", help="Продукт амплификации ОБЯЗАТЕЛЬНО будет содержать этот фрагмент.")
        excluded_reg = rc2.text_input("Исключить регион (Excluded)", placeholder="Например: 200, 30", help="Праймеры НЕ будут ложиться на эту область.")
        
        st.divider()
        st.markdown("**Координаты экзонов (1-based)**")
        exon_text_val = "\n".join([f"{s}-{e}" for s, e in st.session_state.exons_1b]) if st.session_state.exons_1b else ""
        exon_text = st.text_area("Список (один экзон на строку: start-end)", value=exon_text_val, height=100)
        
        if st.button("Сохранить экзоны вручную"):
            try:
                st.session_state.exons_1b = parse_exons_1based(exon_text)
                st.success(f"Экзоны сохранены: {len(st.session_state.exons_1b)}")
            except Exception as e:
                st.error(str(e))

    if seq_input:
        sequences = parse_fasta(seq_input)
        if len(sequences) > 1:
            seq_names = [s[0] for s in sequences]
            selected_name = st.selectbox("Выберите последовательность", seq_names)
            raw_seq = dict(sequences)[selected_name]
            st.session_state.seq_name = selected_name
        elif len(sequences) == 1:
            raw_seq = sequences[0][1]
            if not st.session_state.seq_name or st.session_state.seq_name == "input_sequence":
                 st.session_state.seq_name = sequences[0][0]
        else:
            raw_seq = ""

        is_valid, cleaned_seq, error_msg = validate_sequence(raw_seq)

        if not is_valid:
            st.error(f"❌ {error_msg}")
        else:
            st.session_state.sequence = cleaned_seq
            st.markdown(f"<div class='stat-box'>✅ Готово к подбору. Длина: <b>{len(cleaned_seq)} bp</b> | Матрица: <b>{st.session_state.seq_name}</b> | GC: <b>{calc_gc(cleaned_seq):.1f}%</b></div>", unsafe_allow_html=True)

            if st.button("🚀 Подобрать праймеры", type="primary", use_container_width=True):
                # Парсинг регионов
                t_start, t_len = parse_region(target_reg)
                e_start, e_len = parse_region(excluded_reg)

                params = DesignParams(
                    primer_opt_size=int(primer_opt_size), primer_min_size=int(primer_min_size), primer_max_size=int(primer_max_size),
                    primer_opt_tm=float(primer_opt_tm), primer_min_tm=float(primer_min_tm), primer_max_tm=float(primer_max_tm), primer_max_tm_diff=float(primer_max_tm_diff),
                    primer_min_gc=float(gc_range[0]), primer_max_gc=float(gc_range[1]), primer_gc_clamp=int(primer_gc_clamp),
                    product_size_min=int(product_min), product_size_max=int(product_max), num_return=int(num_return), is_qpcr=(mode == "qPCR по инструкции (RT-qPCR)"),
                    
                    primer_salt_monovalent=float(primer_salt_mono), primer_salt_divalent=float(primer_salt_div),
                    primer_dna_conc=float(primer_dna_conc), primer_dntp_conc=float(primer_dntp_conc),
                    primer_max_self_any=float(max_self_any), primer_max_self_end=float(max_self_end), primer_max_poly_x=int(max_poly_x),
                    
                    target_start=t_start, target_length=t_len, excluded_start=e_start, excluded_length=e_len
                )
                
                with st.spinner("Запуск Primer3..."):
                    errors = get_design_errors(cleaned_seq, params)
                    if errors:
                        st.error(f"❌ {errors}")
                    else:
                        primers = design_primers(cleaned_seq, params)
                        if mode == "qPCR по инструкции (RT-qPCR)" and enforce_junction and st.session_state.exons_1b:
                            junctions = exon_junctions_0based(st.session_state.exons_1b)
                            ann = annotate_primers_with_junctions(primers, junctions, int(min_tail))
                            ann_filtered = [x for x in ann if (x["fwd_spans_junction"] or x["rev_spans_junction"])]
                            
                            filtered_primers = []
                            for p in primers:
                                for a in ann_filtered:
                                    if p.pair_id == a["pair_id"]:
                                        filtered_primers.append(p)
                                        break
                                        
                            st.session_state.primers = filtered_primers
                            
                            df_temp = primers_to_dataframe(filtered_primers)
                            if not df_temp.empty:
                                df_temp["Exon Junction"] = "Да"
                            st.session_state.primers_df = df_temp
                            
                            if ann_filtered: st.success(f"✅ Найдено {len(ann_filtered)} пар с перекрытием exon–exon junction.")
                            else: st.warning("⚠️ Пары с перекрытием стыка экзонов не найдены. Смягчите условия.")
                        else:
                            st.session_state.primers = primers
                            st.session_state.primers_df = primers_to_dataframe(primers)
                            st.success(f"✅ Успешно подобрано {len(primers)} пар праймеров!")

            # ОТОБРАЖЕНИЕ ТАБЛИЦЫ
            if st.session_state.primers_df is not None and not st.session_state.primers_df.empty:
                st.subheader("Результаты подбора (Топ кандидатов)")
                
                df_disp = st.session_state.primers_df.copy()
                
                if "Forward (5'→3')" in df_disp.columns:
                    rename_map = {
                        "Forward (5'→3')": "left_seq", "Reverse (5'→3')": "right_seq", "Продукт (п.н.)": "product_size",
                        "Fwd Tm (°C)": "left_tm", "Rev Tm (°C)": "right_tm", "Fwd GC%": "left_gc", "Rev GC%": "right_gc",
                        "Fwd Self Any": "left_self_any", "Fwd 3' Self": "left_self_end", 
                        "Rev Self Any": "right_self_any", "Rev 3' Self": "right_self_end",
                        "Pair Compl Any": "pair_compl_any", "Pair 3' Compl": "pair_compl_end",
                        "Штраф": "penalty", "#": "pair_id"
                    }
                    df_disp = df_disp.rename(columns=rename_map)
                    df_disp["pair_id"] -= 1

                cols_to_keep = [
                    "pair_id", "left_seq", "right_seq", "product_size", 
                    "left_tm", "right_tm", "left_gc", "right_gc", 
                    "left_self_any", "left_self_end", "right_self_any", "right_self_end",
                    "pair_compl_any", "pair_compl_end", "penalty"
                ]
                if "fwd_spans_junction" in df_disp.columns: cols_to_keep.extend(["fwd_spans_junction", "rev_spans_junction"])

                final_cols = [c for c in cols_to_keep if c in df_disp.columns]
                df_disp = df_disp[final_cols]
                
                final_rename_map = {
                    "pair_id": "Пара #", "left_seq": "Forward", "right_seq": "Reverse", "product_size": "Продукт (bp)", 
                    "left_tm": "F-Tm", "right_tm": "R-Tm", "left_gc": "F-GC%", "right_gc": "R-GC%", "penalty": "Штраф",
                    "left_self_any": "F-Self", "left_self_end": "F-3'Self", 
                    "right_self_any": "R-Self", "right_self_end": "R-3'Self", 
                    "pair_compl_any": "Pair-Compl", "pair_compl_end": "Pair-3'Compl",
                    "fwd_spans_junction": "F-Exon", "rev_spans_junction": "R-Exon"
                }
                df_disp = df_disp.rename(columns=final_rename_map)
                df_disp["Пара #"] += 1
                
                st.dataframe(
                    df_disp,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Штраф": st.column_config.NumberColumn(format="%.2f"),
                        "F-Tm": st.column_config.NumberColumn(format="%.1f °C"),
                        "R-Tm": st.column_config.NumberColumn(format="%.1f °C"),
                        "F-GC%": st.column_config.ProgressColumn(format="%.1f%%", min_value=0, max_value=100),
                        "R-GC%": st.column_config.ProgressColumn(format="%.1f%%", min_value=0, max_value=100),
                        "Продукт (bp)": st.column_config.NumberColumn(format="%d"),
                        "F-Self": st.column_config.NumberColumn(format="%.2f"),
                        "F-3'Self": st.column_config.NumberColumn(format="%.2f"),
                        "R-Self": st.column_config.NumberColumn(format="%.2f"),
                        "R-3'Self": st.column_config.NumberColumn(format="%.2f"),
                        "Pair-Compl": st.column_config.NumberColumn(format="%.2f"),
                        "Pair-3'Compl": st.column_config.NumberColumn(format="%.2f")
                    },
                )


# ─── Вкладка 2: BLAST ──────────────────────────────────────────
with tab_blast:
    st.header("🎯 BLAST & Графическое выравнивание")
    st.markdown("Проверка специфичности праймеров по базам мРНК и Геномной ДНК.")

    blast_source = st.radio("Источник праймеров для проверки:", ["Из списка подобранных", "Ввести вручную (свои праймеры)"], horizontal=True)
    
    blast_left_seq = ""
    blast_right_seq = ""

    if blast_source == "Из списка подобранных":
        if not st.session_state.primers:
            st.info("ℹ️ Сначала подберите праймеры на вкладке «Подбор» или переключитесь на ручной ввод.")
        else:
            pair_options = [f"Пара #{p.pair_id+1} | {p.product_size} п.н. | F: {p.left_seq} | R: {p.right_seq}" for p in st.session_state.primers]
            sel_idx = st.selectbox("Выберите пару для проверки", range(len(pair_options)), format_func=lambda i: pair_options[i])
            selected_pair = st.session_state.primers[sel_idx]

            blast_left_seq = selected_pair.left_seq
            blast_right_seq = selected_pair.right_seq

            col1, col2 = st.columns(2)
            col1.markdown(f'<div class="primer-box"><strong>Forward:</strong><br>{blast_left_seq}</div>', unsafe_allow_html=True)
            col2.markdown(f'<div class="primer-box"><strong>Reverse:</strong><br>{blast_right_seq}</div>', unsafe_allow_html=True)
            
    else: 
        col1, col2 = st.columns(2)
        with col1:
            raw_left = st.text_input("Forward праймер (5'→3')", placeholder="Например: GCGCACAACCGGCAGAA")
            blast_left_seq = re.sub(r"\s+", "", raw_left).upper()
        with col2:
            raw_right = st.text_input("Reverse праймер (5'→3')", placeholder="Например: CCAACCATGACGCCCTGA")
            blast_right_seq = re.sub(r"\s+", "", raw_right).upper()
            
        if blast_left_seq and blast_right_seq:
            st.success("✅ Праймеры введены. Можно запускать проверку.")

    st.divider()

    if blast_left_seq and blast_right_seq:
        col1, col2 = st.columns(2)
        with col1:
            default_dbs = [db for db in ["refseq_rna", "refseq_representative_genomes"] if db in BLAST_DATABASES]
            selected_dbs_keys = st.multiselect("Базы данных", list(BLAST_DATABASES.keys()), default=default_dbs)
            selected_dbs = [BLAST_DATABASES[k] if isinstance(BLAST_DATABASES[k], str) else BLAST_DATABASES[k]["db"] for k in selected_dbs_keys]
        with col2:
            orgs = list(ORGANISM_GROUPS.keys())
            selected_orgs_keys = st.multiselect("Организм", orgs, default=["Человек"] if "Человек" in orgs else [orgs[0]])

        c3, c4 = st.columns(2)
        max_hits = c3.number_input("Макс. хитов", 5, 100, 20)
        identity_threshold = c4.number_input("Порог идентичности (%)", 50.0, 100.0, 80.0, 1.0)

        if st.button("🔍 Запустить BLAST", type="primary", use_container_width=True):
            if not selected_dbs or not selected_orgs_keys:
                st.error("Выберите базу данных и организм.")
            else:
                progress_bar, status_text = st.progress(0), st.empty()
                def blast_progress(p, t): progress_bar.progress(p); status_text.text(t)
                
                with st.spinner("Поиск в NCBI... Это займет время (1-3 мин)"):
                    st.session_state.blast_results = check_primer_pair_specificity(
                        blast_left_seq, blast_right_seq, selected_dbs, selected_orgs_keys, int(max_hits), float(identity_threshold), blast_progress
                    )
                progress_bar.progress(1.0)
                status_text.text("✅ Поиск завершен!")

        if st.session_state.blast_results:
            st.markdown("### 📊 Результаты выравнивания")
            for key, pair_results in st.session_state.blast_results.items():
                st.markdown(f"#### 🗄️ База: `{key}`")
                
                # --- ПОИСК И ВЫВОД ГЕНОВ-МИШЕНЕЙ (АМПЛИКОНОВ) ---
                if len(pair_results) == 2:
                    fwd_res, rev_res = pair_results[0], pair_results[1]
                    
                    # Группируем хиты по Accession ID для поиска пересечений
                    from collections import defaultdict
                    fwd_dict = defaultdict(list)
                    for h in fwd_res.hits: fwd_dict[h.accession].append(h)
                    
                    rev_dict = defaultdict(list)
                    for h in rev_res.hits: rev_dict[h.accession].append(h)
                    
                    # Пересечение: на этот accession сели ОБА праймера
                    common_accs = set(fwd_dict.keys()).intersection(set(rev_dict.keys()))
                    
                    if common_accs:
                        st.markdown("##### 🎯 Возможные гены-мишени (Ампликоны)")
                        for acc in common_accs:
                            f_hit = fwd_dict[acc][0]  # Берем лучший хит
                            r_hit = rev_dict[acc][0]
                            
                            # Вычисляем теоретический размер ПЦР-продукта
                            coords = [f_hit.subject_start, f_hit.subject_end, r_hit.subject_start, r_hit.subject_end]
                            prod_size = max(coords) - min(coords) + 1
                            
                            st.success(
                                f"**Ген:** {f_hit.gene_description}\n\n"
                                f"🧬 **Accession:** `{acc}` | **Размер продукта:** `~{prod_size} bp` | **Организм:** *{f_hit.organism}*"
                            )
                    else:
                        st.warning("⚠️ Праймеры не образуют совместного продукта (ампликона) на известных матрицах в этой базе.")
                        
                st.divider()

                # Вывод детального выравнивания для каждого праймера
                for res in pair_results:
                    if res.error_message:
                        st.error(f"**{res.primer_name}**: {res.error_message}")
                        continue
                    if res.warning_message: 
                        st.warning(f"**{res.primer_name}**: {res.warning_message}")
                    else: 
                        st.success(f"**{res.primer_name}**: Высокоспецифичен.")

                    for h in res.hits[:10]:
                        with st.expander(f"➤ {h.accession} | {h.organism} | Идент: {h.identity}% | E-val: {h.e_value:.2e}"):
                            st.caption(f"Описание: {h.title}")
                            html_alignment = render_blast_alignment_html(
                                getattr(h, "query_seq", ""), 
                                getattr(h, "subject_seq", ""), 
                                h.query_start, 
                                h.subject_start, 
                                getattr(h, "match_seq", None)
                            )
                            st.markdown(html_alignment, unsafe_allow_html=True)


# ─── Вкладка 3: In-Silico PCR ─────────────────────────────────
with tab_manual:
    st.header("🛠️ In-Silico PCR & Калькулятор праймеров")
    
    st.markdown("### 🧮 1. Ручной анализ праймеров (Термодинамика)")
    col1, col2 = st.columns(2)
    manual_fwd = col1.text_input("Forward праймер (для анализа)", placeholder="ATGCGTACG...")
    manual_rev = col2.text_input("Reverse праймер (для анализа)", placeholder="CGTACG...")

    if st.button("Рассчитать параметры", use_container_width=True):
        valid_fwd, clean_fwd, _ = validate_sequence(manual_fwd) if manual_fwd else (False, "", "")
        valid_rev, clean_rev, _ = validate_sequence(manual_rev) if manual_rev else (False, "", "")

        c1, c2, c3, c4 = st.columns(4)
        if valid_fwd:
            c1.metric("F-Tm", f"{calc_tm(clean_fwd):.1f} °C")
            c2.metric("F-GC", f"{calc_gc(clean_fwd):.1f} %")
            c3.metric("F-Hairpin", f"{calc_hairpin_tm(clean_fwd):.1f} °C")
            c4.metric("F-Dimer", f"{calc_homodimer_tm(clean_fwd):.1f} °C")
        if valid_rev:
            c1.metric("R-Tm", f"{calc_tm(clean_rev):.1f} °C")
            c2.metric("R-GC", f"{calc_gc(clean_rev):.1f} %")
            c3.metric("R-Hairpin", f"{calc_hairpin_tm(clean_rev):.1f} °C")
            c4.metric("R-Dimer", f"{calc_homodimer_tm(clean_rev):.1f} °C")
        if valid_fwd and valid_rev:
            het_tm = calc_heterodimer_tm(clean_fwd, clean_rev)
            st.metric("Гетеродимер F+R", f"{het_tm:.1f} °C", delta="Критично!" if het_tm > min(calc_tm(clean_fwd), calc_tm(clean_rev))-10 else "Ок", delta_color="inverse")

    st.divider()

    st.markdown("### 🔬 2. Виртуальная ПЦР (In-Silico PCR)")
    st.markdown("Ищет возможные продукты на заданной матрице, допуская указанное количество несовпадений (мисматчей).")
    
    vpcr_template = st.text_area("Матрица для виртуальной ПЦР (Sequence)", value=st.session_state.sequence if st.session_state.sequence else "", height=150)
    col_v1, col_v2 = st.columns(2)
    with col_v1:
        vpcr_max_mm = st.number_input("Макс. мисматчей (на каждый праймер)", 0, 5, 2)
    with col_v2:
        vpcr_max_size = st.number_input("Макс. размер продукта (п.н.)", 50, 20000, 5000)

    if st.button("🧬 Запустить виртуальную ПЦР", type="primary", use_container_width=True):
        if not vpcr_template or not manual_fwd or not manual_rev:
            st.error("Укажите матрицу и оба праймера (в блоке 1 выше)!")
        else:
            _, t_clean, _ = validate_sequence(vpcr_template)
            _, f_clean, _ = validate_sequence(manual_fwd)
            _, r_clean, _ = validate_sequence(manual_rev)
            
            with st.spinner("Симуляция ПЦР..."):
                products = run_virtual_pcr(t_clean, f_clean, r_clean, max_mismatches=int(vpcr_max_mm), max_product=int(vpcr_max_size))
            
            if not products:
                st.warning("⚠️ Продукты амплификации не найдены. Попробуйте увеличить допуск по мисматчам.")
            else:
                st.success(f"✅ Найдено возможных продуктов: {len(products)}")
                for i, prod in enumerate(products, 1):
                    is_danger = prod['fwd_3prime_mm'] > 0 or prod['rev_3prime_mm'] > 0
                    bg_color = "#FFF3E0" if is_danger else "#E8F5E9"
                    border_color = "#FF9800" if is_danger else "#4CAF50"
                    
                    st.markdown(f"""
                    <div style="background:{bg_color}; border-left:5px solid {border_color}; padding:15px; border-radius:8px; margin-bottom:10px; color:#202124;">
                        <h4 style="margin-top:0; color:#1a73e8;">Продукт #{i} — <b>{prod['size']} bp</b></h4>
                        <ul style="margin-bottom:0;">
                            <li><b>Локация:</b> {prod['fwd_start']} ... {prod['rev_end']}</li>
                            <li><b>Forward мисматчи:</b> {prod['fwd_mismatches']} (из них на 3'-конце: <span style="color:{'red' if prod['fwd_3prime_mm']>0 else 'green'}"><b>{prod['fwd_3prime_mm']}</b></span>)</li>
                            <li><b>Reverse мисматчи:</b> {prod['rev_mismatches']} (из них на 3'-конце: <span style="color:{'red' if prod['rev_3prime_mm']>0 else 'green'}"><b>{prod['rev_3prime_mm']}</b></span>)</li>
                        </ul>
                        {f"<p style='color:#d32f2f; margin-bottom:0; margin-top:10px; font-weight:bold;'>⚠️ Внимание: Есть мисматчи на 3'-концах! Вероятность амплификации сильно снижена.</p>" if is_danger else ""}
                    </div>
                    """, unsafe_allow_html=True)


# ─── Вкладка 4: Визуализация и экспорт ────────────────────────
with tab_results:
    st.header("📊 Визуализация и Экспорт")
    if not st.session_state.primers:
        st.info("Сначала подберите праймеры на вкладке «Подбор».")
    else:
        primers = st.session_state.primers
        seq_length = len(st.session_state.sequence) if st.session_state.sequence else 1000

        pairs_to_show = st.multiselect("Выберите пары на карте", list(range(len(primers))), default=list(range(min(5, len(primers)))), format_func=lambda i: f"Пара #{i+1}")
        
        if pairs_to_show:
            fig_map = create_primer_map(seq_length, primers, pairs_to_show)
            if st.session_state.exons_1b: 
                fig_map = add_exons_to_map(fig_map, st.session_state.exons_1b)
            st.plotly_chart(fig_map, use_container_width=True)

        col1, col2 = st.columns(2)
        with col1: st.plotly_chart(create_tm_comparison_chart(primers), use_container_width=True)
        with col2: st.plotly_chart(create_gc_comparison_chart(primers), use_container_width=True)
        
        fig_pen = create_penalty_chart(primers)
        st.plotly_chart(fig_pen, use_container_width=True)

        st.divider()
        st.subheader("📥 Экспорт результатов")

        col1, col2, col3 = st.columns(3)

        with col1:
            if st.session_state.primers_df is not None:
                st.download_button(
                    "⬇️ Скачать таблицу (CSV)", 
                    data=export_to_csv(st.session_state.primers_df), 
                    file_name=f"primers_{st.session_state.seq_name}.csv", 
                    mime="text/csv", 
                    use_container_width=True
                )
        with col2:
            mode_str = "qpcr" if mode == "qPCR по инструкции (RT-qPCR)" else "standard"
            report_text = generate_primer_report(primers, sequence_name=st.session_state.seq_name, mode=mode_str)
            st.download_button(
                "⬇️ Скачать текстовый отчет (TXT)", 
                data=report_text, 
                file_name=f"primer_report_{st.session_state.seq_name}.txt", 
                mime="text/plain",
                use_container_width=True
            )
        with col3:
            order_lines = ["Name\tSequence\tScale\tPurification"]
            for p in primers:
                order_lines.append(f"Fwd_{p.pair_id+1}\t{p.left_seq}\t25nm\tSTD")
                order_lines.append(f"Rev_{p.pair_id+1}\t{p.right_seq}\t25nm\tSTD")
            order_data = "\n".join(order_lines)
            st.download_button(
                "⬇️ Файл для заказа (TSV)",
                data=order_data,
                file_name=f"order_{st.session_state.seq_name}.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
            )

# ─── Футер ──────────────────────────────────────────────────────
st.divider()
st.markdown("""
<div style="text-align:center; color:#999; font-size:0.85rem;">
    🧬 PrimerDesigner v2.0 | Полная интеграция NCBI, Графический BLAST и In-Silico PCR | 
    <a href="https://primer3.org" target="_blank" style="color:#1a73e8;">Primer3</a> · 
    <a href="https://blast.ncbi.nlm.nih.gov" target="_blank" style="color:#1a73e8;">NCBI BLAST</a>
</div>
""", unsafe_allow_html=True)