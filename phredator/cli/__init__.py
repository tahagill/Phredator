"""CLI package

This package exposes a `main` symbol so entry points that refer to
`phredator.cli:main` continue to work when the CLI implementation lives
in `phredator.cli.cli`.
"""

from .cli import main

__all__ = ["main"]
