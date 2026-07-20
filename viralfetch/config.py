"""Configuration resolution: NCBI email / API key, cache & config paths.

Precedence for a value is: explicit CLI flag > environment variable >
persisted config file. There is deliberately **no default email** — the NCBI
policy requires a real one, and inventing a value is forbidden (SPEC section 3).
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

from platformdirs import user_cache_dir, user_config_dir

APP_NAME = "viralfetch"

CONFIG_DIR = Path(user_config_dir(APP_NAME))
CONFIG_FILE = CONFIG_DIR / "config.json"
CACHE_DIR = Path(user_cache_dir(APP_NAME))


class ConfigError(Exception):
    """Raised for missing-but-required configuration (e.g. no NCBI email)."""


@dataclass
class Config:
    email: str | None = None
    api_key: str | None = None
    format: str = "rich"  # "rich" | "json"
    verbose: bool = False
    no_cache: bool = False

    def require_email(self) -> str:
        """Return the NCBI email or fail loudly — never fabricate one."""
        if not self.email:
            raise ConfigError(
                "No NCBI email configured. Set $NCBI_EMAIL, pass --email, or run "
                "`viralfetch config --store-ncbi-email you@example.com`."
            )
        return self.email

    @property
    def rate_limit(self) -> int:
        """Allowed NCBI requests/second (3 without an API key, 10 with one)."""
        return 10 if self.api_key else 3


def _load_file() -> dict:
    if CONFIG_FILE.is_file():
        try:
            return json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def resolve(
    *,
    email: str | None = None,
    api_key: str | None = None,
    fmt: str = "rich",
    verbose: bool = False,
    no_cache: bool = False,
) -> Config:
    """Resolve effective configuration from flags, environment, and file."""
    stored = _load_file()
    return Config(
        email=email or os.environ.get("NCBI_EMAIL") or stored.get("email"),
        api_key=api_key or os.environ.get("NCBI_API_KEY") or stored.get("api_key"),
        format=fmt,
        verbose=verbose,
        no_cache=no_cache,
    )


def store(*, email: str | None = None, api_key: str | None = None) -> Path:
    """Persist email and/or API key to the config file. Returns its path."""
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    data = _load_file()
    if email is not None:
        data["email"] = email
    if api_key is not None:
        data["api_key"] = api_key
    CONFIG_FILE.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return CONFIG_FILE
