"""
Модуль для работы с NCBI API (Entrez).
Позволяет скачивать GenBank-файлы по Accession ID и извлекать
нуклеотидную последовательность вместе с координатами экзонов.
"""
from Bio import Entrez, SeqIO
from typing import Tuple, List, Optional

# Укажите email (требование NCBI для использования их API)
Entrez.email = "primer_designer_app@example.com"

def fetch_sequence_and_exons(accession: str) -> Tuple[Optional[str], List[Tuple[int, int]], str]:
    """
    Скачивает GenBank запись по Accession ID.
    
    Returns:
        (sequence, exons_1based, error_message)
    """
    accession = accession.strip()
    if not accession:
        return None, [], "Accession ID не может быть пустым."

    try:
        # efetch загружает данные из базы nucleotide
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        sequence = str(record.seq).upper()
        exons_1based = []

        # Парсим аннотации для поиска экзонов
        for feature in record.features:
            if feature.type == "exon":
                # В Biopython индексы 0-based и полуоткрытые [start:end)
                # Переводим в 1-based включительные для нашего интерфейса
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                exons_1based.append((start, end))

        # Сортируем экзоны по стартовой позиции
        exons_1based.sort(key=lambda x: x[0])

        return sequence, exons_1based, ""

    except Exception as e:
        return None, [], f"Ошибка при загрузке из NCBI (возможно, неверный ID): {str(e)}"
