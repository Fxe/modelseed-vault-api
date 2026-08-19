from typing import Any
import stat
from pathlib import Path
import os
import tomllib
import atexit


class ConfigError(RuntimeError):
    """Raised for anything wrong with the config file or a profile."""


class BioSEED:

    def __init__(self, default_config=None):
        self.default_config = default_config
        if default_config is None:
            self.default_config = Path(
                os.environ.get("BIOSEED_CONFIG", "~/.config/bioseed/config.toml")
            ).expanduser()

    def profiles(self) -> dict[str, list[str]]:
        """What's configured, by client kind. Handy for tab-completion-less discovery."""
        return {kind: sorted(v) for kind, v in self._config().items() if isinstance(v, dict)}

    def _config(self) -> dict[str, Any]:
        path = self.default_config
        if not path.exists():
            raise ConfigError(
                f"No config at {path}. Run `python -m bioseed.db init` to create one."
            )
        mode = path.stat().st_mode
        if mode & (stat.S_IRWXG | stat.S_IRWXO):
            raise ConfigError(f"{path} is group- or world-accessible. Run: chmod 600 {path}")
        with path.open("rb") as fh:
            return tomllib.load(fh)

    def _profile(self, kind: str, name: str) -> dict[str, Any]:
        section = self._config().get(kind, {})
        if name not in section:
            available = ", ".join(sorted(section)) or "none configured"
            raise ConfigError(
                f"No [{kind}.{name}] in {self.default_config}. Available {kind} profiles: {available}"
            )
        return dict(section[name])

    def _secret(self, kind: str, name: str, entry: dict[str, Any]) -> str:
        """Resolve a password. Order: env override, password_env, password_command, literal."""
        override = f"BIOSEED_{kind}_{name}_PASSWORD".upper().replace("-", "_")
        if override in os.environ:
            return os.environ[override]

        if "password_env" in entry:
            var = entry["password_env"]
            print(var)
            if var not in os.environ:
                raise ConfigError(
                    f"[{kind}.{name}] expects the password in ${var}, which is unset. "
                    f"Set it, or set ${override} instead."
                )
            return os.environ[var]

        if "password_command" in entry:
            cmd = entry["password_command"]
            try:
                out = subprocess.run(
                    shlex.split(cmd), capture_output=True, text=True, check=True, timeout=30
                )
            except subprocess.CalledProcessError as exc:
                raise ConfigError(
                    f"password_command for [{kind}.{name}] failed: {exc.stderr.strip()}"
                ) from None
            return out.stdout.strip()

        if "password" in entry:
            return entry["password"]

        raise ConfigError(
            f"[{kind}.{name}] has no password source. Add password_env, "
            f"password_command, or password."
        )

    def mysql(self, name: str = "default", **engine_kwargs: Any):
        """Return a cached SQLAlchemy Engine for the named MySQL profile."""
        from sqlalchemy import create_engine
        from sqlalchemy.engine import URL

        p = self._profile("mysql", name)
        url = URL.create(
            "mysql+pymysql",
            username=p["user"],
            password=self._secret("mysql", name, p),
            host=p["host"],
            port=p.get("port", 3306),
            database=p.get("database", name),
        )
        opts: dict[str, Any] = {
            "pool_pre_ping": True,  # recover from a connection the server already dropped
            "pool_recycle": 3600,  # stay under MySQL wait_timeout on long-lived kernels
            "pool_size": 5,
            "max_overflow": 5,
        }
        opts.update(engine_kwargs)
        engine = create_engine(url, **opts)
        atexit.register(engine.dispose)
        return engine

    def mongo(self, name: str = "poplar"):
        from pymongo import MongoClient

        p = self._profile("mongo", name)
        mongo_client = MongoClient(p["host"], int(p["port"]))
        return mongo_client[p["database"]]

    def vault(self, name: str = "poplar"):
        from modelseed_vault.vault import Vault
        p = self._profile("vault", name)
        vault = Vault(p["host"])
        return vault

    def minio(self, name: str = "poplar"):
        from minio import Minio

        p = self._profile("minio", name)
        host = p['host']
        port = p['port']
        client = Minio(
            f"{host}:{port}",
            access_key=p['user'],
            secret_key=self._secret("minio", name, p),
            secure=False,  # Set to True if using HTTPS
        )

        # atexit.register()

        return client
