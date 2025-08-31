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

    # Skip the first 9 tokens (equivalent to 6 + 3 extractions in C++)
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
            x = float(a)  # noqa: F841 (x is unused; kept for parity with C++)
            y = float(b)
        except ValueError:
            # C++ >> extraction would fail here and break out of the loop
            break
        pairs.append((x, y))

    if not pairs:
        sys.stderr.write(f"{missing_label} has no numeric data after header\n")
        sys.exit(1)

    # Return the last y-value (equivalent to disX[naX-1][1] in the C++ code)
    return pairs[-1][1]


def main() -> None:
    # Map input files to their "inputFileN" labels for C++-style error messages
    y2 = read_last_y("output1.txt", "inputFile1")  # corresponds to dis2 in C++
    y1 = read_last_y("output2.txt", "inputFile2")  # corresponds to dis1 in C++
    y3 = read_last_y("output3.txt", "inputFile3")
    y4 = read_last_y("output4.txt", "inputFile4")

    # Mimic default C++ iostream formatting (precision = 6, general format)
    def fmt(v: float) -> str:
        return f"{v:.6g}"

    line = f"{fmt(y2)}    {fmt(y1)}    {fmt(y3)}    {fmt(y4)}\n"
    with open("output_all.txt", "w") as out:
        out.write(line)


if __name__ == "__main__":
    main()
