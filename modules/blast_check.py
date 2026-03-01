"""
Модуль проверки специфичности праймеров через NCBI BLAST API.
Реализован кастомный алгоритм поллинга (Polling) для обхода тайм-аутов
облачных серверов (Streamlit Cloud).
"""
import time
import requests
from io import StringIO
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple
from Bio.Blast import NCBIXML

# Полный список групп организмов
ORGANISM_GROUPS: Dict[str, Dict[str, str]] = {
    "Все организмы": {"taxid": "", "label": "Все организмы (без фильтра)"},
    "Человек": {"taxid": "9606", "label": "Homo sapiens (TaxID: 9606)"},
    "Мышь": {"taxid": "10090", "label": "Mus musculus (TaxID: 10090)"},
    "Крыса": {"taxid": "10116", "label": "Rattus norvegicus (TaxID: 10116)"},
    "Бактерии": {"taxid": "2", "label": "Bacteria (TaxID: 2)"},
    "Археи": {"taxid": "2157", "label": "Archaea (TaxID: 2157)"},
    "Эукариоты": {"taxid": "2759", "label": "Eukaryota (TaxID: 2759)"},
    "Простейшие": {"taxid": "2759", "label": "Protozoa (подмножество Eukaryota)"},
    "Грибы": {"taxid": "4751", "label": "Fungi (TaxID: 4751)"},
    "Растения": {"taxid": "3193", "label": "Embryophyta / Viridiplantae (TaxID: 3193)"},
    "Насекомые": {"taxid": "50557", "label": "Insecta (TaxID: 50557)"},
    "Вирусы": {"taxid": "10239", "label": "Viruses (TaxID: 10239)"},
}

BLAST_DATABASES: Dict[str, str] = {
    "refseq_rna": "RefSeq mRNA (для проверки целевого транскрипта)",
    "refseq_representative_genomes": "RefSeq Genomes (для проверки на геномную ДНК)",
    "nt": "Nucleotide collection (nt) — общая база",
    "nr": "Non-redundant protein (nr) — для blastp",
}

@dataclass
class BlastHit:
    accession: str
    title: str
    organism: str
    identity: float 
    coverage: float 
    e_value: float
    alignment_length: int
    mismatches: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    score: float
    query_seq: str = ""
    match_seq: str = ""
    subject_seq: str = ""
    gene_description: str = ""

@dataclass
class SpecificityResult:
    primer_name: str
    primer_seq: str
    database: str
    organism_filter: str
    total_hits: int = 0
    hits: List[BlastHit] = field(default_factory=list)
    is_specific: bool = True
    warning_message: str = ""
    error_message: str = ""
    target_genes: set = field(default_factory=set)


def run_blast_check(
    primer_seq: str,
    primer_name: str = "primer",
    database: str = "nt",
    organism_group: str = "Все организмы",
    max_hits: int = 20,
    evalue_threshold: float = 1000.0,
    identity_threshold: float = 80.0,
    max_retries: int = 3,
    progress_callback=None,
    fraction: float = 0.0,
    total_text: str = ""
) -> SpecificityResult:
    
    result = SpecificityResult(
        primer_name=primer_name,
        primer_seq=primer_seq,
        database=database,
        organism_filter=organism_group,
    )

    taxid = ORGANISM_GROUPS.get(organism_group, {}).get("taxid", "")
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

    put_params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": database,
        "QUERY": primer_seq,
        "SHORT_QUERY": "true",
        "EXPECT": str(evalue_threshold),
        "HITLIST_SIZE": str(max_hits),
        "FORMAT_TYPE": "XML",
        "EMAIL": "primer_designer_pro@science.com",
        "TOOL": "PrimerDesignerPro",
        "WORD_SIZE": "7"
    }
    if taxid:
        put_params["ENTREZ_QUERY"] = f"txid{taxid}[ORGN]"

    for attempt in range(max_retries):
        try:
            if progress_callback:
                atmpt_txt = f" (Попытка {attempt+1})" if attempt > 0 else ""
                progress_callback(fraction, f"{total_text}{atmpt_txt} - Отправка в NCBI...")

            # 1. Отправляем запрос (Put)
            resp = requests.post(url, data=put_params, timeout=30)
            resp.raise_for_status()

            rid = None
            for line in resp.text.splitlines():
                if line.startswith("    RID = "):
                    rid = line.split("=")[1].strip()
                    break

            if not rid:
                raise Exception("Не удалось получить RID (уникальный ID запроса) от NCBI.")

            # 2. Ждём готовности (Get/SearchInfo) - ИМЕННО ЭТО СПАСАЕТ ОТ ЗАВИСАНИЙ
            poll_params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
            start_time = time.time()
            xml_data = None

            while True:
                time.sleep(5)
                elapsed = int(time.time() - start_time)
                if progress_callback:
                    # Постоянное обновление UI не дает облаку сбросить соединение
                    progress_callback(fraction, f"{total_text} - NCBI ищет (RID: {rid})... {elapsed} сек.")

                poll_resp = requests.get(url, params=poll_params, timeout=30)
                poll_text = poll_resp.text
                
                if "Status=WAITING" in poll_text:
                    continue
                if "Status=FAILED" in poll_text:
                    raise Exception("Поиск NCBI завершился с ошибкой (FAILED).")
                if "Status=UNKNOWN" in poll_text:
                    raise Exception("Неизвестный статус RID. Сессия устарела.")
                if "Status=READY" in poll_text:
                    if "There Are No Hits" in poll_text:
                        xml_data = "" # Нет хитов
                    break

            if xml_data == "": 
                return result # Возвращаем 0 хитов без ошибок

            # 3. Скачиваем результат
            if progress_callback:
                progress_callback(fraction, f"{total_text} - Скачивание результатов (RID: {rid})...")

            ret_params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
            ret_resp = requests.get(url, params=ret_params, timeout=60)
            ret_resp.raise_for_status()
            xml_data = ret_resp.text

            # 4. Парсим XML
            blast_records = NCBIXML.parse(StringIO(xml_data))
            try:
                blast_record = next(blast_records)
            except StopIteration:
                return result

            for alignment in blast_record.alignments:
                title = alignment.title or "Unknown Title"
                clean_title = title.split('|')[-1].strip() if '|' in title else title
                if clean_title.startswith("PREDICTED: "):
                    clean_title = clean_title.replace("PREDICTED: ", "")
                
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100
                    coverage_pct = (abs(hsp.query_end - hsp.query_start) + 1) / len(primer_seq) * 100

                    organism = ""
                    if "[" in title and "]" in title:
                        organism = title[title.rfind("[")+1:title.rfind("]")]

                    hit = BlastHit(
                        accession=alignment.accession,
                        title=title[:200],
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
                        query_seq=hsp.query,    
                        match_seq=hsp.match,    
                        subject_seq=hsp.sbjct,
                        gene_description=clean_title
                    )
                    result.hits.append(hit)
                    if identity_pct >= identity_threshold and coverage_pct >= 80:
                        result.target_genes.add(clean_title)

            result.total_hits = len(result.hits)

            high_identity_hits = [h for h in result.hits if h.identity >= identity_threshold and h.coverage >= 80]
            if len(high_identity_hits) > 5:
                result.is_specific = False
                result.warning_message = f"Найдено {len(high_identity_hits)} хитов с идентичностью ≥{identity_threshold}% и покрытием ≥80%. Праймер неспецифичен!"
            elif len(high_identity_hits) > 1:
                result.warning_message = "Найдено несколько совпадений. Проверьте список генов-мишеней (вероятно, это сплайс-варианты)."

            return result

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 429:
                err_str = "NCBI временно заблокировал доступ (Слишком много запросов). Попробуйте позже."
            else:
                err_str = f"HTTP Ошибка {e.response.status_code}: {str(e)}"
        except Exception as e:
            err_str = str(e)
            if not err_str: err_str = repr(e)
                
        if attempt < max_retries - 1:
            time.sleep(3)
            continue
        else:
            result.error_message = f"Сбой API NCBI: {err_str}"
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
    
    results = {}
    total_checks = len(databases) * len(organism_groups) * 2
    current = 0

    for db in databases:
        for org in organism_groups:
            key = f"{db} | {org}"
            pair_results = []

            for primer_seq, primer_name in [(left_seq, "Forward"), (right_seq, "Reverse")]:
                current += 1
                fraction = current / total_checks
                text = f"BLAST: {primer_name} → {db} [{org}] ({current}/{total_checks})"
                
                if progress_callback:
                    progress_callback(fraction, text)

                spec_result = run_blast_check(
                    primer_seq=primer_seq,
                    primer_name=primer_name,
                    database=db,
                    organism_group=org,
                    max_hits=max_hits,
                    identity_threshold=identity_threshold,
                    progress_callback=progress_callback,
                    fraction=fraction,
                    total_text=text
                )
                pair_results.append(spec_result)

            results[key] = pair_results

    return results

def format_blast_summary(results: Dict[str, List[SpecificityResult]]) -> str:
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