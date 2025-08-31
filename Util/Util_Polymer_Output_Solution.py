#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

def read_last_y(filename: str, missing_label: str) -> float:
    try:
        with open(filename, "r") as f:
            tokens = f.read().split()
    except OSError:
        sys.stderr.write(f"{missing_label} could not be opened\n")
        sys.exit(1)

    # Skip first 9 tokens (header)
    data_tokens = tokens[9:] if len(tokens) > 9 else []

    # Parse numeric pairs until failure
    pairs = []
    it = iter(data_tokens)
    while True:
        try:
            a = next(it)
            b = next(it)
        except StopIteration:
            break
        try:
            x = float(a)  # kept for parity with C++ but unused
            y = float(b)
        except ValueError:
            # Emulate stream extraction failure (stop reading pairs)
            break
        pairs.append((x, y))

    # In the original C++, if there are no valid pairs this would lead to
    # undefined indexing. Here we raise a clear error to avoid silent failure.
    if not pairs:
        raise RuntimeError(f"{missing_label}: no numeric (x, y) pairs after header")

    return pairs[-1][1]


def fmt(v: float) -> str:
    # Mimic default C++ iostream "general" with precision ~6
    return f"{v:.6g}"


def main() -> None:
    y2 = read_last_y("output1.txt", "inputFile1")  # dis2 from output1.txt
    y1 = read_last_y("output2.txt", "inputFile2")  # dis1 from output2.txt
    y3 = read_last_y("output3.txt", "inputFile3")
    y4 = read_last_y("output4.txt", "inputFile4")

    line = f"{fmt(y2)}    {fmt(y1)}    {fmt(y3)}    {fmt(y4)}    \n"
    with open("output_all.txt", "w") as out:
        out.write(line)


if __name__ == "__main__":
    main()
