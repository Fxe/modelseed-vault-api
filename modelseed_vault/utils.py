import hashlib


def sha_hex(s: str):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()
