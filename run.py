#!/usr/bin/env python3
"""
PrimerPro - Приложение для подбора праймеров (PCR / qPCR)
с проверкой специфичности по базам данных NCBI BLAST.

Запуск: python run.py
"""
import subprocess
import sys
import os

def check_dependencies():
    """Проверка и установка зависимостей."""
    required = {
        "streamlit": "streamlit",
        "primer3": "primer3-py",
        "Bio": "biopython",
        "requests": "requests",
        "pandas": "pandas",
        "plotly": "plotly",
    }
    missing = []
    for module, package in required.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)

    if missing:
        print(f"[PrimerPro] Устанавливаю зависимости: {', '.join(missing)}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", *missing, "-q"]
        )
        print("[PrimerPro] Зависимости установлены!")

def main():
    check_dependencies()
    app_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py")
    print("[PrimerPro] Запускаю приложение...")
    subprocess.run([sys.executable, "-m", "streamlit", "run", app_path,
                    "--server.headless", "true"])

if __name__ == "__main__":
    main()