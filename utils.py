"""
Утилиты: валидация последовательностей, экспорт результатов,
расчёт термодинамических свойств, In-Silico PCR.
"""
import re
import primer3
import pandas as pd
from typing import List, Optional, Tuple, Dict
from modules.primer_design import PrimerResult


def validate_sequence(seq: str) -> Tuple[bool, str, str]:
    """Проверяет корректность нуклеотидной последовательности."""
    cleaned = re.sub(r"[\s\d>]", "", seq.upper())
    if cleaned.startswith(">"):
        lines = seq.strip().split("\n")
        cleaned = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
        cleaned = re.sub(r"[\s\d]", "", cleaned.upper())

    if not cleaned:
        return False, "", "Последовательность пуста."
    if len(cleaned) < 50:
        return False, cleaned, "Последовательность слишком короткая (минимум 50 п.н.)."
    if len(cleaned) > 200000:
        return False, cleaned, "Последовательность слишком длинная (максимум 200 000 п.н.)."

    invalid = set(cleaned) - set("ATGCNRYSWKMBDHV")
    if invalid:
        return False, cleaned, f"Недопустимые символы: {', '.join(invalid)}"

    n_fraction = cleaned.count("N") / len(cleaned)
    if n_fraction > 0.1:
        return False, cleaned, f"Слишком много неопределённых нуклеотидов (N): {n_fraction*100:.1f}%"

    return True, cleaned, ""


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    sequences = []
    current_name = ""
    current_seq = []

    for line in text.strip().split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_name and current_seq:
                sequences.append((current_name, "".join(current_seq)))
            current_name = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(re.sub(r"[\s\d]", "", line.upper()))

    if current_name and current_seq:
        sequences.append((current_name, "".join(current_seq)))

    if not sequences:
        cleaned = re.sub(r"[\s\d]", "", text.upper())
        if cleaned:
            sequences.append(("input_sequence", cleaned))

    return sequences


def calc_tm(seq: str) -> float:
    try: return round(primer3.calc_tm(seq), 1)
    except: return 0.0

def calc_gc(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 1) if seq else 0.0

def calc_hairpin_tm(seq: str) -> float:
    try: return round(primer3.calc_hairpin(seq).tm, 1)
    except: return 0.0

def calc_homodimer_tm(seq: str) -> float:
    try: return round(primer3.calc_homodimer(seq).tm, 1)
    except: return 0.0

def calc_heterodimer_tm(seq1: str, seq2: str) -> float:
    try: return round(primer3.calc_heterodimer(seq1, seq2).tm, 1)
    except: return 0.0

def reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
                  "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M",
                  "M": "K", "B": "V", "V": "B", "D": "H", "H": "D"}
    return "".join(complement.get(b, "N") for b in reversed(seq.upper()))


def primers_to_dataframe(primers: List[PrimerResult]) -> pd.DataFrame:
    """Конвертирует список PrimerResult в DataFrame с учетом комплементарности."""
    rows = []
    for p in primers:
        row = {
            "#": p.pair_id + 1,
            "Forward (5'→3')": p.left_seq,
            "Reverse (5'→3')": p.right_seq,
            "Fwd Tm (°C)": round(p.left_tm, 1),
            "Rev Tm (°C)": round(p.right_tm, 1),
            "Fwd GC%": round(p.left_gc, 1),
            "Rev GC%": round(p.right_gc, 1),
            "Fwd длина": p.left_length,
            "Rev длина": p.right_length,
            "Продукт (п.н.)": p.product_size,
            
            # Добавленные параметры комплементарности для таблицы
            "Fwd Self Any": round(p.left_self_any, 2) if p.left_self_any is not None else 0.0,
            "Fwd 3' Self": round(p.left_self_end, 2) if p.left_self_end is not None else 0.0,
            "Rev Self Any": round(p.right_self_any, 2) if p.right_self_any is not None else 0.0,
            "Rev 3' Self": round(p.right_self_end, 2) if p.right_self_end is not None else 0.0,
            "Pair Compl Any": round(p.pair_compl_any, 2) if p.pair_compl_any is not None else 0.0,
            "Pair 3' Compl": round(p.pair_compl_end, 2) if p.pair_compl_end is not None else 0.0,
            
            "Штраф": round(p.penalty, 2),
            "Fwd позиция": p.left_start,
            "Rev позиция": p.right_start,
        }
        if p.probe_seq:
            row["Проба"] = p.probe_seq
            row["Проба Tm"] = round(p.probe_tm, 1) if p.probe_tm else ""
            row["Проба GC%"] = round(p.probe_gc, 1) if p.probe_gc else ""
            
        rows.append(row)

    return pd.DataFrame(rows)


def export_to_csv(df: pd.DataFrame) -> str:
    return df.to_csv(index=False)

def generate_primer_report(primers: List[PrimerResult], sequence_name: str = "input", mode: str = "standard") -> str:
    lines = [
        f"# Отчёт о подборе праймеров",
        f"Последовательность: {sequence_name}",
        f"Режим: {'qPCR' if mode == 'qpcr' else 'Стандартный ПЦР'}",
        f"Найдено пар: {len(primers)}",
        "=" * 60,
    ]
    for p in primers:
        lines.append(f"\n## Пара #{p.pair_id + 1} (штраф: {p.penalty:.2f})")
        lines.append(f"  Forward:  5\'-{p.left_seq}-3\'  (Tm={p.left_tm:.1f}°C, GC={p.left_gc:.1f}%, {p.left_length} п.н.)")
        lines.append(f"  Reverse:  5\'-{p.right_seq}-3\'  (Tm={p.right_tm:.1f}°C, GC={p.right_gc:.1f}%, {p.right_length} п.н.)")
        lines.append(f"  Продукт: {p.product_size} п.н.")
        
        if p.pair_compl_any is not None:
            lines.append(f"  Комплементарность пары (any): {p.pair_compl_any:.1f}")
            lines.append(f"  Комплементарность пары (end): {p.pair_compl_end:.1f}")
            
    return "\n".join(lines)


# ─── IN-SILICO PCR (Виртуальная ПЦР) ───────────────────────────

def find_binding_sites(template: str, search_seq: str, max_mismatches: int = 2, is_reverse: bool = False) -> List[Dict]:
    """Ищет места посадки скользящим окном с учетом мисматчей."""
    sites = []
    t_len = len(template)
    p_len = len(search_seq)
    if t_len == 0 or p_len == 0 or t_len < p_len:
        return sites

    search_upper = search_seq.upper()
    template_upper = template.upper()

    for i in range(t_len - p_len + 1):
        window = template_upper[i:i+p_len]
        
        # Быстрый подсчет несовпадений
        mismatches = sum(1 for a, b in zip(window, search_upper) if a != b)
        
        if mismatches <= max_mismatches:
            # Считаем мисматчи на 3'-конце праймера (последние 5 нуклеотидов).
            # Для Forward праймера это КОНЕЦ поисковой строки.
            # Для Reverse праймера (ищем обратно-комплементарную) 3'-конец праймера - это НАЧАЛО поисковой строки!
            if not is_reverse:
                end_mismatches = sum(1 for a, b in zip(window[-5:], search_upper[-5:]) if a != b)
            else:
                end_mismatches = sum(1 for a, b in zip(window[:5], search_upper[:5]) if a != b)

            sites.append({
                "start": i + 1,  # 1-based
                "end": i + p_len,
                "mismatches": mismatches,
                "3prime_mismatches": end_mismatches,
                "sequence": window
            })
    return sites

def run_virtual_pcr(template: str, fwd: str, rev: str, max_mismatches: int = 2, max_product: int = 5000) -> List[Dict]:
    """Запускает симуляцию ПЦР, находя продукты амплификации."""
    template = template.upper()
    rev_comp = reverse_complement(rev)

    fwd_sites = find_binding_sites(template, fwd, max_mismatches, is_reverse=False)
    rev_sites = find_binding_sites(template, rev_comp, max_mismatches, is_reverse=True)

    products = []
    for f in fwd_sites:
        for r in rev_sites:
            # Праймер Reverse должен быть ПОСЛЕ праймера Forward
            if f["start"] < r["end"]:
                product_size = r["end"] - f["start"] + 1
                if product_size <= max_product:
                    products.append({
                        "fwd_start": f["start"],
                        "fwd_end": f["end"],
                        "rev_start": r["start"],
                        "rev_end": r["end"],
                        "size": product_size,
                        "fwd_mismatches": f["mismatches"],
                        "rev_mismatches": r["mismatches"],
                        "fwd_3prime_mm": f["3prime_mismatches"],
                        "rev_3prime_mm": r["3prime_mismatches"]
                    })
    
    # Сортируем по размеру продукта
    products.sort(key=lambda x: x["size"])
    return products