"""Progress logging without import-time filesystem side effects."""
from datetime import datetime
from pathlib import Path


def print_to_log(logfile, message, warning=False):
    level = "WARNING | " if warning else ""
    text = f"{datetime.now():%Y-%m-%d %H:%M:%S} | {level}{message}"
    with Path(logfile).open("a", encoding="utf-8") as stream:
        stream.write(text + "\n")
    print(text, flush=True)
