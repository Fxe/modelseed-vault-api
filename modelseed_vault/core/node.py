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
