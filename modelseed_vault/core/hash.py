import hashlib


def _hash_string(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


class HashString(str):

    def __new__(cls, v, strip_ending_star=True):
        # print('validate!!', v)
        if strip_ending_star:
            v = HashString.strip_ending_star(v)
        instance = super().__new__(cls, v.upper())
        return instance

    @property
    def hash_value(self):
        h = _hash_string(self)
        return h

    @staticmethod
    def strip_ending_star(s):
        if s.endswith('*'):
            return s[:-1]
        return s


class HashStringList(list):

    def append(self, o, /):
        if type(o) is str:
            super().append(HashString(o))
        elif type(o) is HashString:
            super().append(o)
        else:
            raise ValueError('bad type')

    @property
    def hash_value(self):
        h_list = [x.hash_value for x in self]
        hash_seq = "_".join(sorted(h_list))
        return _hash_string(hash_seq)
