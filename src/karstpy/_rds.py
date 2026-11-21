"""Write the R data types used by the ALFA-K export on Python 3.10."""

from dataclasses import dataclass
import gzip
import struct

import numpy as np
import pandas as pd


@dataclass
class RAttributed:
    value: object
    attributes: dict


class _Writer:
    NIL, SYMBOL, PAIR, CHAR = 254, 1, 2, 9
    LOGICAL, INTEGER, REAL, STRINGS, VECTOR = 10, 13, 14, 16, 19
    OBJECT, ATTRIBUTES, TAG = 1 << 8, 1 << 9, 1 << 10

    def __init__(self, stream):
        self.stream = stream

    def integers(self, *values):
        self.stream.write(struct.pack(f">{len(values)}i", *values))

    def string(self, value):
        if value is None or pd.isna(value):
            self.integers(self.CHAR, -1)
            return
        encoded = str(value).encode("utf-8")
        ascii_flag, utf8_flag = 64, 8
        flags = ascii_flag if encoded.isascii() else utf8_flag
        self.integers(self.CHAR | (flags << 12), len(encoded))
        self.stream.write(encoded)

    def attributes(self, attributes):
        for name, value in attributes.items():
            self.integers(self.PAIR | self.TAG, self.SYMBOL)
            self.string(name)
            self.write(value)
        self.integers(self.NIL)

    def write(self, value, attributes=None):
        attributes = {} if attributes is None else attributes.copy()
        if isinstance(value, RAttributed):
            self.write(value.value, value.attributes)
            return
        if value is None:
            self.integers(self.NIL)
            return
        if isinstance(value, pd.DataFrame):
            attributes.update(
                {"class": "data.frame", "row.names": np.array([-(2**31), -len(value)], dtype=np.int32)
                 if len(value) else np.empty(0, dtype=np.int32)}
            )
            self.write({name: column.to_numpy() for name, column in value.items()}, attributes)
            return
        if isinstance(value, dict):
            attributes["names"] = np.array(list(value), dtype=str)
            value = list(value.values())
        flags = self.ATTRIBUTES if attributes else 0
        if "class" in attributes:
            flags |= self.OBJECT
        if isinstance(value, (list, tuple)):
            self.integers(self.VECTOR | flags, len(value))
            for item in value:
                self.write(item)
        else:
            array = np.asarray(value)
            if array.ndim > 1:
                attributes["dim"] = np.array(array.shape, dtype=np.int32)
                flags |= self.ATTRIBUTES
            array = array.ravel(order="F")
            kind = array.dtype.kind
            if kind in "OUS":
                self.integers(self.STRINGS | flags, len(array))
                for item in array:
                    self.string(item)
            elif kind in "bui":
                self.integers((self.LOGICAL if kind == "b" else self.INTEGER) | flags, len(array))
                limits = np.iinfo(np.int32)
                if np.any(array < limits.min) or np.any(array > limits.max):
                    raise OverflowError("R integer value is outside the 32-bit range.")
                self.stream.write(array.astype(">i4").tobytes())
            elif kind == "f":
                self.integers(self.REAL | flags, len(array))
                encoded = np.ascontiguousarray(array, dtype=">f8")
                if np.isnan(encoded).any():
                    encoded = encoded.copy()
                    encoded.view(">u8")[np.isnan(encoded)] = 0x7FF00000000007A2
                self.stream.write(encoded.tobytes())
            else:
                raise TypeError(f"Unsupported RDS value with dtype {array.dtype}.")
        if attributes:
            self.attributes(attributes)


def write_rds(path, value):
    """Write a gzip-compressed version 2 XDR RDS stream."""
    with gzip.open(path, "wb") as stream:
        stream.write(b"X\n")
        writer = _Writer(stream)
        writer.integers(2, 0x40402, 0x20300)
        writer.write(value)
