"""
Модуль подбора праймеров с использованием primer3-py.
Поддерживает обычные праймеры, строго настроенные qPCR-специфичные параметры,
а также расширенные настройки химии (соли, концентрации) и регионов.
"""
import primer3
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class PrimerResult:
    """Результат подбора одной пары праймеров."""
    pair_id: int
    left_seq: str
    right_seq: str
    left_tm: float
    right_tm: float
    left_gc: float
    right_gc: float
    left_start: int
    left_length: int
    right_start: int
    right_length: int
    product_size: int
    penalty: float
    # qPCR-специфичные
    product_tm: Optional[float] = None
    left_self_any: Optional[float] = None
    left_self_end: Optional[float] = None
    right_self_any: Optional[float] = None
    right_self_end: Optional[float] = None
    pair_compl_any: Optional[float] = None
    pair_compl_end: Optional[float] = None
    # Опциональный probe для qPCR (TaqMan)
    probe_seq: Optional[str] = None
    probe_tm: Optional[float] = None
    probe_gc: Optional[float] = None


@dataclass
class DesignParams:
    """Параметры для дизайна праймеров."""
    # Основные
    primer_opt_size: int = 20
    primer_min_size: int = 18
    primer_max_size: int = 27
    primer_opt_tm: float = 60.0
    primer_min_tm: float = 57.0
    primer_max_tm: float = 63.0
    primer_max_tm_diff: float = 3.0
    primer_min_gc: float = 40.0
    primer_max_gc: float = 60.0
    primer_opt_gc_percent: float = 50.0
    primer_gc_clamp: int = 0
    
    # Продукт
    product_size_min: int = 100
    product_size_max: int = 1000
    
    # Концентрации и соли (влияют на расчет Tm)
    primer_salt_monovalent: float = 50.0  # mM (K+, Na+)
    primer_salt_divalent: float = 1.5     # mM (Mg2+)
    primer_dna_conc: float = 50.0         # nM (концентрация праймера)
    primer_dntp_conc: float = 0.6         # mM (dNTPs)
    
    # Вторичные структуры и повторы
    primer_max_poly_x: int = 4
    primer_max_self_any: float = 8.0
    primer_max_self_end: float = 3.0
    primer_max_ns_accepted: int = 0
    
    # Количество
    num_return: int = 10
    
    # qPCR-режим
    is_qpcr: bool = False
    pick_probe: bool = False  # для TaqMan
    overlap_junction: Optional[int] = None # Позиция стыка экзонов
    
    # Целевой регион (опционально)
    target_start: Optional[int] = None
    target_length: Optional[int] = None
    
    # Исключаемый регион
    excluded_start: Optional[int] = None
    excluded_length: Optional[int] = None


def get_qpcr_defaults() -> DesignParams:
    """Возвращает параметры, строго оптимизированные для RT-qPCR по протоколу."""
    return DesignParams(
        primer_opt_size=20,
        primer_min_size=15,
        primer_max_size=25,
        primer_opt_tm=60.0,
        primer_min_tm=58.0,
        primer_max_tm=62.0,
        primer_max_tm_diff=2.0,   # Max Tm diff = 2-3
        primer_min_gc=45.0,       # GC 45-60%
        primer_max_gc=60.0,
        primer_opt_gc_percent=52.5,
        primer_gc_clamp=2,        # GC-clamp = 2
        product_size_min=100,     # Product size 100-300
        product_size_max=300,
        primer_salt_monovalent=50.0,
        primer_salt_divalent=3.0, # В qPCR-смесях Mg2+ часто повышен
        primer_dna_conc=250.0,    # Концентрация праймеров в qPCR обычно 200-400 nM
        primer_dntp_conc=0.8,     # В смесях обычно по 0.2 mM каждого (итого 0.8)
        primer_max_poly_x=3,      # Строже к повторам в qPCR
        primer_max_self_any=6.0,  # Max Self Any = 6.00
        primer_max_self_end=2.0,  # Max 3' Self = 2.00
        num_return=10,
        is_qpcr=True,
        pick_probe=False,
    )


def get_standard_defaults() -> DesignParams:
    """Возвращает стандартные параметры для обычного ПЦР."""
    return DesignParams()


def design_primers(sequence: str, params: DesignParams) -> List[PrimerResult]:
    """
    Подбирает праймеры для заданной последовательности.
    """
    # Очистка последовательности
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    # Формируем seq_args
    seq_args: Dict[str, Any] = {
        "SEQUENCE_ID": "input_sequence",
        "SEQUENCE_TEMPLATE": sequence,
    }

    if params.target_start is not None and params.target_length is not None:
        seq_args["SEQUENCE_TARGET"] = [params.target_start, params.target_length]

    if params.excluded_start is not None and params.excluded_length is not None:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [[params.excluded_start, params.excluded_length]]
        
    # Принудительное перекрытие стыка экзонов
    if params.overlap_junction is not None:
        seq_args["SEQUENCE_OVERLAP_JUNCTION_LIST"] = [params.overlap_junction]

    # Формируем global_args
    global_args: Dict[str, Any] = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 1 if params.pick_probe else 0,
        "PRIMER_NUM_RETURN": params.num_return,
        "PRIMER_OPT_SIZE": params.primer_opt_size,
        "PRIMER_MIN_SIZE": params.primer_min_size,
        "PRIMER_MAX_SIZE": params.primer_max_size,
        "PRIMER_OPT_TM": params.primer_opt_tm,
        "PRIMER_MIN_TM": params.primer_min_tm,
        "PRIMER_MAX_TM": params.primer_max_tm,
        "PRIMER_PAIR_MAX_DIFF_TM": params.primer_max_tm_diff,
        "PRIMER_MIN_GC": params.primer_min_gc,
        "PRIMER_MAX_GC": params.primer_max_gc,
        "PRIMER_OPT_GC_PERCENT": params.primer_opt_gc_percent,
        "PRIMER_GC_CLAMP": params.primer_gc_clamp,
        "PRIMER_PRODUCT_SIZE_RANGE": [[params.product_size_min, params.product_size_max]],
        
        # Химия и соли
        "PRIMER_SALT_MONOVALENT": params.primer_salt_monovalent,
        "PRIMER_SALT_DIVALENT": params.primer_salt_divalent,
        "PRIMER_DNA_CONC": params.primer_dna_conc,
        "PRIMER_DNTP_CONC": params.primer_dntp_conc,
        
        # Настройки самокомплементарности и повторов
        "PRIMER_MAX_SELF_ANY": params.primer_max_self_any,
        "PRIMER_MAX_SELF_END": params.primer_max_self_end,
        "PRIMER_MAX_POLY_X": params.primer_max_poly_x,
        "PRIMER_MAX_NS_ACCEPTED": params.primer_max_ns_accepted,
        
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": 1,
    }

    if params.is_qpcr:
        global_args["PRIMER_PRODUCT_OPT_SIZE"] = (params.product_size_min + params.product_size_max) // 2

    # Запуск Primer3
    result = primer3.design_primers(seq_args, global_args)

    # Парсинг результатов
    primers = []
    num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)

    for i in range(num_returned):
        left_pos = result.get(f"PRIMER_LEFT_{i}", (0, 0))
        right_pos = result.get(f"PRIMER_RIGHT_{i}", (0, 0))

        pr = PrimerResult(
            pair_id=i,
            left_seq=result.get(f"PRIMER_LEFT_{i}_SEQUENCE", ""),
            right_seq=result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", ""),
            left_tm=result.get(f"PRIMER_LEFT_{i}_TM", 0.0),
            right_tm=result.get(f"PRIMER_RIGHT_{i}_TM", 0.0),
            left_gc=result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0.0),
            right_gc=result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0.0),
            left_start=left_pos[0],
            left_length=left_pos[1],
            right_start=right_pos[0],
            right_length=right_pos[1],
            product_size=result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0),
            penalty=result.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0),
            product_tm=result.get(f"PRIMER_PAIR_{i}_PRODUCT_TM", None),
            left_self_any=result.get(f"PRIMER_LEFT_{i}_SELF_ANY", None),
            left_self_end=result.get(f"PRIMER_LEFT_{i}_SELF_END", None),
            right_self_any=result.get(f"PRIMER_RIGHT_{i}_SELF_ANY", None),
            right_self_end=result.get(f"PRIMER_RIGHT_{i}_SELF_END", None),
            # ДОБАВЛЕНО ИЗВЛЕЧЕНИЕ КОМПЛЕМЕНТАРНОСТИ ПАРЫ
            pair_compl_any=result.get(f"PRIMER_PAIR_{i}_COMPL_ANY", None),
            pair_compl_end=result.get(f"PRIMER_PAIR_{i}_COMPL_END", None),
        )

        if params.pick_probe:
            pr.probe_seq = result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", None)
            pr.probe_tm = result.get(f"PRIMER_INTERNAL_{i}_TM", None)
            pr.probe_gc = result.get(f"PRIMER_INTERNAL_{i}_GC_PERCENT", None)

        primers.append(pr)

    return primers


def get_design_errors(sequence: str, params: DesignParams) -> Optional[str]:
    """Возвращает ошибки дизайна, если они есть."""
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    seq_args = {
        "SEQUENCE_ID": "input_sequence",
        "SEQUENCE_TEMPLATE": sequence,
    }
    if params.target_start is not None and params.target_length is not None:
        seq_args["SEQUENCE_TARGET"] = [params.target_start, params.target_length]
        
    if params.excluded_start is not None and params.excluded_length is not None:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [[params.excluded_start, params.excluded_length]]
        
    if params.overlap_junction is not None:
        seq_args["SEQUENCE_OVERLAP_JUNCTION_LIST"] = [params.overlap_junction]

    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_NUM_RETURN": 1,
        "PRIMER_OPT_SIZE": params.primer_opt_size,
        "PRIMER_MIN_SIZE": params.primer_min_size,
        "PRIMER_MAX_SIZE": params.primer_max_size,
        "PRIMER_OPT_TM": params.primer_opt_tm,
        "PRIMER_MIN_TM": params.primer_min_tm,
        "PRIMER_MAX_TM": params.primer_max_tm,
        "PRIMER_PAIR_MAX_DIFF_TM": params.primer_max_tm_diff,
        "PRIMER_MIN_GC": params.primer_min_gc,
        "PRIMER_MAX_GC": params.primer_max_gc,
        "PRIMER_GC_CLAMP": params.primer_gc_clamp,
        "PRIMER_PRODUCT_SIZE_RANGE": [[params.product_size_min, params.product_size_max]],
        "PRIMER_SALT_MONOVALENT": params.primer_salt_monovalent,
        "PRIMER_SALT_DIVALENT": params.primer_salt_divalent,
        "PRIMER_DNA_CONC": params.primer_dna_conc,
        "PRIMER_DNTP_CONC": params.primer_dntp_conc,
        "PRIMER_MAX_SELF_ANY": params.primer_max_self_any,
        "PRIMER_MAX_SELF_END": params.primer_max_self_end,
        "PRIMER_MAX_POLY_X": params.primer_max_poly_x,
        "PRIMER_MAX_NS_ACCEPTED": params.primer_max_ns_accepted,
    }

    result = primer3.design_primers(seq_args, global_args)

    errors = result.get("PRIMER_ERROR", "")
    warnings = result.get("PRIMER_WARNING", "")

    if errors:
        return f"Ошибка: {errors}"
    if result.get("PRIMER_PAIR_NUM_RETURNED", 0) == 0:
        explain = result.get("PRIMER_PAIR_EXPLAIN", "Нет подходящих пар")
        left_explain = result.get("PRIMER_LEFT_EXPLAIN", "")
        right_explain = result.get("PRIMER_RIGHT_EXPLAIN", "")
        msg = f"Не найдено подходящих праймеров.\nПары: {explain}"
        if left_explain:
            msg += f"\nLeft: {left_explain}"
        if right_explain:
            msg += f"\nRight: {right_explain}"
        if warnings:
            msg += f"\nПредупреждения: {warnings}"
        return msg

    return None