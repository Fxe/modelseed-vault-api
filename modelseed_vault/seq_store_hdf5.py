import hashlib
import h5py


class SeqStoreHdf5:
    def __init__(self, h5_store):
        self.CHARSET_VALIDATION = {}
        self.h5_store = h5_store
        self.string_dtype = h5py.string_dtype(encoding="utf-8")
        pass

    def validate_sequence(self, s):
        v = set(s) - self.CHARSET_VALIDATION
        if len(v) == 0:
            return True
        raise ValueError("Invalid input - unaccepted characters:" + ", ".join(v))

    @staticmethod
    def get_sequence_hash(seq):
        return hashlib.sha256(seq.encode("utf-8")).hexdigest()

    def store_sequence(self, seq):
        self.validate_sequence(seq)
        h = self.get_sequence_hash(seq)
        if h not in self.h5_store:
            ds = self.h5_store.create_dataset(
                h, (1,), dtype=self.string_dtype, compression="gzip"
            )
            ds[0] = seq
        return h

    def get_sequence(self, hash_val):
        if hash_val in self.h5_store:
            return self.h5_store[hash_val][0].decode("utf-8")
        return None

    def close(self):
        self.h5_store.close()
        return self.h5_store


def load_dna_seq_store_h5(path, mode):
    h5_store = h5py.File(path, mode)
    s = SeqStoreHdf5(h5_store)
    # K	Guanine / Thymine
    # M	Adenine / Cytosine
    # S	Guanine / Cytosine
    # W	Adenine / Thymine
    # B	Guanine / Thymine / Cytosine
    # D	Guanine / Adenine / Thymine
    # H	Adenine / Cytosine / Thymine
    # V	Guanine / Cytosine / Adenine
    s.CHARSET_VALIDATION = {
        "A",
        "C",
        "G",
        "T",
        "U",  # Adenine, Cytosine, Guanine, Thymine, Uracil
        "R",
        "Y",  # Guanine / Adenine (purine), Cytosine / Thymine (pyrimidine)
        "K",
        "M",
        "S",
        "W",
        "B",
        "D",
        "H",
        "V",
        "N",  # Adenine / Guanine / Cytosine / Thymine
    }
    return s


def load_protein_seq_store_h5(path, mode):
    h5_store = h5py.File(path, mode)
    s = SeqStoreHdf5(h5_store)
    s.CHARSET_VALIDATION = {
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    }
    return s
