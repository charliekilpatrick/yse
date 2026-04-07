"""Shared logging setup for CLIs and scripts.

Do **not** name this module ``logging`` — that would shadow :mod:`logging` and break
``import logging`` in unexpected ways. Use :func:`configure_logging` once at program
entry (e.g. in ``main()``).

**Streams**

* **stdout** — :const:`logging.DEBUG` and :const:`logging.INFO` (progress and detail).
* **stderr** — :const:`logging.WARNING` and above (issues that should not be mistaken
  for normal progress output).
"""

from __future__ import annotations

import logging
import sys
from typing import Final

__all__ = ["configure_logging", "DEFAULT_FORMAT"]

DEFAULT_FORMAT: Final[str] = "%(asctime)s %(levelname)s [%(name)s] %(message)s"
_DEFAULT_DATEFMT: Final[str] = "%H:%M:%S"


class _MaxLevelFilter(logging.Filter):
    """Pass records with ``levelno <= max_level`` (keeps DEBUG/INFO off stderr)."""

    def __init__(self, max_level: int) -> None:
        super().__init__()
        self.max_level = max_level

    def filter(self, record: logging.LogRecord) -> bool:
        return record.levelno <= self.max_level


def configure_logging(
    *,
    level: int = logging.INFO,
    fmt: str = DEFAULT_FORMAT,
    datefmt: str = _DEFAULT_DATEFMT,
) -> None:
    """Attach stdout (INFO and below) and stderr (WARNING+) handlers to the root logger.

    Clears existing root handlers so repeated CLI invocations in one process do not
    duplicate output. Idempotent for typical single-shot ``main()`` use.

    Parameters
    ----------
    level
        Minimum level emitted (e.g. :const:`logging.DEBUG` with ``--verbose``).
    fmt, datefmt
        Passed to :class:`logging.Formatter`.
    """
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(level)

    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setLevel(logging.DEBUG)
    stdout_h.addFilter(_MaxLevelFilter(logging.INFO))
    stdout_h.setFormatter(logging.Formatter(fmt, datefmt=datefmt))

    stderr_h = logging.StreamHandler(sys.stderr)
    stderr_h.setLevel(logging.WARNING)
    stderr_h.setFormatter(logging.Formatter(fmt, datefmt=datefmt))

    root.addHandler(stdout_h)
    root.addHandler(stderr_h)

    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
