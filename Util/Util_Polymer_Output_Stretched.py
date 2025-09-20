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

    if len(tokens) <= 9:
        sys.stderr.write(f"{missing_label} has insufficient tokens after header\n")
        sys.exit(1)
    data_tokens = tokens[9:]

    # Read numeric pairs (x, y) until encountering non-numeric content
    pairs = []
    it = iter(data_tokens)
    while True:
        try:
            a = next(it)
            b = next(it)
        except StopIteration:
            break
        try:
            x = float(a)
            y = float(b)
        except ValueError:
            break
        pairs.append((x, y))

    if not pairs:
        sys.stderr.write(f"{missing_label} has no numeric data after header\n")
        sys.exit(1)

    return pairs[-1][1]


def main() -> None:
    y2 = read_last_y("output1.txt", "inputFile1")
    y1 = read_last_y("output2.txt", "inputFile2")
    y3 = read_last_y("output3.txt", "inputFile3")
    y4 = read_last_y("output4.txt", "inputFile4")

    def fmt(v: float) -> str:
        return f"{v:.6g}"

    line = f"{fmt(y2)}    {fmt(y1)}    {fmt(y3)}    {fmt(y4)}\n"
    with open("output_all.txt", "w") as out:
        out.write(line)


if __name__ == "__main__":
    main()
