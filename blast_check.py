"""
Модуль проверки специфичности праймеров через NCBI BLAST API.
Адаптировано для корректных проверок RT-qPCR (мРНК + геномка).
"""
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple
from io import StringIO

from Bio.Blast import NCBIWWW, NCBIXML


# Предустановленные группы организмов с NCBI TaxID
ORGANISM_GROUPS: Dict[str, Dict[str, str]] = {
    "Человек": {"taxid": "9606", "label": "Homo sapiens (TaxID: 9606)"},
    "Мышь": {"taxid": "10090", "label": "Mus musculus (TaxID: 10090)"},
    "Крыса": {"taxid": "10116", "label": "Rattus norvegicus (TaxID: 10116)"},
    "Бактерии": {"taxid": "2", "label": "Bacteria (TaxID: 2)"},
    "Эукариоты": {"taxid": "2759", "label": "Eukaryota (TaxID: 2759)"},
    "Все организмы": {"taxid": "", "label": "Все организмы (без фильтра)"},
}

# Доступные BLAST-базы (приоритетные для qPCR вынесены наверх)
BLAST_DATABASES: Dict[str, str] = {
    "refseq_rna": "RefSeq mRNA (для проверки целевого транскрипта)",
    "refseq_representative_genomes": "RefSeq Genomes (для проверки на геномную ДНК)",
    "nt": "Nucleotide collection (nt) — общая база",
    "nr": "Non-redundant protein (nr) — для blastp",
}


@dataclass
class BlastHit:
    """Один хит BLAST для праймера."""
    accession: str
    title: str
    organism: str
    identity: float  # процент идентичности
    coverage: float  # покрытие запроса
    e_value: float
    alignment_length: int
    mismatches: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    score: float
    # Поля для графического выравнивания (Добавлено!)
    query_seq: str = ""
    match_seq: str = ""
    subject_seq: str = ""


@dataclass
class SpecificityResult:
    """Результат проверки специфичности одного праймера."""
    primer_name: str
    primer_seq: str
    database: str
    organism_filter: str
    total_hits: int = 0
    hits: List[BlastHit] = field(default_factory=list)
    is_specific: bool = True
    warning_message: str = ""
    error_message: str = ""


def run_blast_check(
    primer_seq: str,
    primer_name: str = "primer",
    database: str = "nt",
    organism_group: str = "Все организмы",
    max_hits: int = 20,
    evalue_threshold: float = 1000,
    identity_threshold: float = 80.0,
    max_retries: int = 3,
) -> SpecificityResult:
    """
    Запускает BLAST-поиск для проверки специфичности праймера.
    """
    result = SpecificityResult(
        primer_name=primer_name,
        primer_seq=primer_seq,
        database=database,
        organism_filter=organism_group,
    )

    taxid = ORGANISM_GROUPS.get(organism_group, {}).get("taxid", "")

    for attempt in range(max_retries):
        try:
            entrez_query = ""
            if taxid:
                entrez_query = f"txid{taxid}[ORGN]"

            blast_handle = NCBIWWW.qblast(
                program="blastn",
                database=database,
                sequence=primer_seq,
                word_size=7,
                expect=evalue_threshold,
                hitlist_size=max_hits,
                entrez_query=entrez_query if entrez_query else None,
                format_type="XML",
            )

            blast_records = NCBIXML.parse(blast_handle)
            blast_record = next(blast_records)

            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100
                    coverage_pct = (abs(hsp.query_end - hsp.query_start) + 1) / len(primer_seq) * 100

                    title = alignment.title
                    organism = ""
                    if "[" in title and "]" in title:
                        organism = title[title.rfind("[")+1:title.rfind("]")]

                    hit = BlastHit(
                        accession=alignment.accession,
                        title=alignment.title[:200],
                        organism=organism,
                        identity=round(identity_pct, 1),
                        coverage=round(coverage_pct, 1),
                        e_value=hsp.expect,
                        alignment_length=hsp.align_length,
                        mismatches=hsp.align_length - hsp.identities,
                        query_start=hsp.query_start,
                        query_end=hsp.query_end,
                        subject_start=hsp.sbjct_start,
                        subject_end=hsp.sbjct_end,
                        score=hsp.score,
                        query_seq=hsp.query,    # Сохраняем выравнивание
                        match_seq=hsp.match,    # Сохраняем выравнивание
                        subject_seq=hsp.sbjct   # Сохраняем выравнивание
                    )
                    result.hits.append(hit)

            result.total_hits = len(result.hits)

            # Оценка специфичности: ищем хиты с высокой гомологией
            high_identity_hits = [h for h in result.hits if h.identity >= identity_threshold and h.coverage >= 80]
            if len(high_identity_hits) > 5:
                result.is_specific = False
                result.warning_message = (
                    f"Найдено {len(high_identity_hits)} хитов с идентичностью ≥{identity_threshold}% "
                    f"и покрытием ≥80%. Праймер может быть неспецифичен!"
                )
            elif len(high_identity_hits) > 1:
                result.warning_message = (
                    f"Найдено {len(high_identity_hits)} совпадений. Убедитесь, что это целевой продукт или его сплайс-варианты."
                )

            return result

        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(5 * (attempt + 1))
                continue
            else:
                result.error_message = f"Ошибка BLAST: {str(e)}"
                return result

    return result


def check_primer_pair_specificity(
    left_seq: str,
    right_seq: str,
    databases: List[str],
    organism_groups: List[str],
    max_hits: int = 20,
    identity_threshold: float = 80.0,
    progress_callback=None,
) -> Dict[str, List[SpecificityResult]]:
    """Проверяет пару праймеров по базам и группам организмов."""
    results = {}
    total_checks = len(databases) * len(organism_groups) * 2
    current = 0

    for db in databases:
        for org in organism_groups:
            key = f"{db} | {org}"
            pair_results = []

            for primer_seq, primer_name in [(left_seq, "Forward"), (right_seq, "Reverse")]:
                current += 1
                if progress_callback:
                    progress_callback(current / total_checks, 
                                     f"BLAST: {primer_name} → {db} [{org}] ({current}/{total_checks})")

                spec_result = run_blast_check(
                    primer_seq=primer_seq,
                    primer_name=primer_name,
                    database=db,
                    organism_group=org,
                    max_hits=max_hits,
                    identity_threshold=identity_threshold,
                )
                pair_results.append(spec_result)
                time.sleep(3)

            results[key] = pair_results

    return results


def format_blast_summary(results: Dict[str, List[SpecificityResult]]) -> str:
    """Форматирует сводку BLAST-результатов."""
    lines = []
    for key, pair in results.items():
        lines.append(f"\n### {key}")
        for res in pair:
            status = "✅ Специфичен" if res.is_specific else "⚠️ Неспецифичен"
            if res.error_message:
                status = f"❌ Ошибка: {res.error_message}"
            lines.append(f"  **{res.primer_name}** ({res.primer_seq}): {status} | Хитов: {res.total_hits}")
            if res.warning_message:
                lines.append(f"    ⚠ {res.warning_message}")
    return "\n".join(lines)