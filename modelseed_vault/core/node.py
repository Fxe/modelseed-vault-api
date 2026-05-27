class Node:
    def __init__(self, key: str, label: str, data=None):
        if not key:
            raise ValueError('empty key')
        if not label:
            raise ValueError('empty label')
        self._key = key.strip()
        self.label = label.strip()
        self.data = data if data else {}

    @property
    def key(self):
        return self._key.replace(" ", "_")

    @property
    def id(self):
        return f"{self.label}/{self.key}"

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return isinstance(other, Node) and self.id == other.id

    def to_json(self):
        out = {
            '_key': self._key
        }
        out.update(self.data)

        return out


class HashNode(Node):

    @staticmethod
    def hash_value(value: str) -> str:
        import hashlib
        return hashlib.sha256(value.encode()).hexdigest()

    def __init__(self, key: str, label: str, data=None):
        hashed_key = self.hash_value(key)
        _data = {}
        if data:
            if '_value' in data:
                if data['_value'] != key:
                    raise ValueError(f"error data contains _value={data['_value']}")
            _data.update(data)

        _data['_value'] = key

        super().__init__(hashed_key, label, data=_data)

