"""
Конфигурация приложения PrimerPro.
"""

# ─── Параметры Primer3 для обычной ПЦР ───
PCR_DEFAULTS = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MAX_SELF_ANY_TH": 45.0,
    "PRIMER_MAX_SELF_END_TH": 35.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0,
    "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[200, 1000]],
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 3.0,
    "PRIMER_DNA_CONC": 250.0,
    "PRIMER_TM_FORMULA": 1,  # SantaLucia 1998
}

# ─── Параметры Primer3 для qPCR ───
QPCR_DEFAULTS = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_MAX_SELF_ANY_TH": 40.0,
    "PRIMER_MAX_SELF_END_TH": 30.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 40.0,
    "PRIMER_PAIR_MAX_COMPL_END_TH": 30.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 200]],
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 3.0,
    "PRIMER_DNA_CONC": 250.0,
    "PRIMER_TM_FORMULA": 1,
    # qPCR-специфичные
    "PRIMER_MAX_HAIRPIN_TH": 47.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
}

# ─── Базы данных NCBI для проверки специфичности ───
BLAST_DATABASES = {
    "nt (все нуклеотиды)":        {"db": "nt",        "description": "Полная нуклеотидная БД NCBI"},
    "refseq_rna (RefSeq RNA)":    {"db": "refseq_rna","description": "Референсные последовательности РНК"},
    "refseq_representative_genomes": {"db": "refseq_representative_genomes",
                                       "description": "Репрезентативные геномы RefSeq"},
}

# Организмы для Entrez-фильтрации
ORGANISM_FILTERS = {
    "Без фильтра":          "",
    "Человек (Homo sapiens)":         "Homo sapiens[ORGN]",
    "Мышь (Mus musculus)":            "Mus musculus[ORGN]",
    "Крыса (Rattus norvegicus)":      "Rattus norvegicus[ORGN]",
    "Бактерии (Bacteria)":            "Bacteria[ORGN]",
    "Простейшие (Protozoa)":          "Protozoa[ORGN]",
    "Грибы (Fungi)":                  "Fungi[ORGN]",
    "Растения (Viridiplantae)":       "Viridiplantae[ORGN]",
    "Насекомые (Insecta)":            "Insecta[ORGN]",
    "Рыбы (Actinopterygii)":         "Actinopterygii[ORGN]",
    "Вирусы (Viruses)":              "Viruses[ORGN]",
    "Архебактерии (Archaea)":         "Archaea[ORGN]",
    "Эукариоты (Eukaryota)":         "Eukaryota[ORGN]",
}

# Цветовая схема
COLORS = {
    "primary": "#0068C9",
    "success": "#21C354",
    "warning": "#FACA2B",
    "error":   "#FF4B4B",
}
