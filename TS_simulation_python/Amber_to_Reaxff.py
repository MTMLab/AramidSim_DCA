#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

INPUT_FILE = "system_after_npt.data"
OUTPUT_FILE = "system.data_modified"

def _read_all_tokens_after_skipping_first_line(fp):
    lines = fp.readlines()
    if not lines:
        return []
    rest = "".join(lines[1:])
    return rest.split()

def _next(tokens, cast=str):
    try:
        t = tokens.pop(0)
    except IndexError:
        raise RuntimeError("Unexpected end of tokens while parsing.")
    if cast is str:
        return t
    try:
        return cast(t)
    except Exception:
        raise RuntimeError(f"Failed to parse token '{t}' as {cast.__name__}")

def main():
    try:
        with open(INPUT_FILE, "r", encoding="utf-8", errors="replace") as f:
            tokens = _read_all_tokens_after_skipping_first_line(f)
    except FileNotFoundError:
        print("Unable to open system.data", file=sys.stderr)
        sys.exit(1)

    kkk = _next(tokens, int); _ = _next(tokens)                       # "<atoms> atoms"
    n_type0 = _next(tokens, int); _ = _next(tokens); _ = _next(tokens) # "<types> atom types"
    n_bond0 = _next(tokens, int); _ = _next(tokens)                    # "<bonds> bonds"
    n_bond_type0 = _next(tokens, int); _ = _next(tokens); _ = _next(tokens)
    n_angle0 = _next(tokens, int); _ = _next(tokens)
    n_angle_type0 = _next(tokens, int); _ = _next(tokens); _ = _next(tokens)
    n_dihedral0 = _next(tokens, int); _ = _next(tokens)
    n_dihedral_type0 = _next(tokens, int); _ = _next(tokens); _ = _next(tokens)
    n_improper0 = _next(tokens, int); _ = _next(tokens)
    n_improper_type0 = _next(tokens, int); _ = _next(tokens); _ = _next(tokens)

    # Box bounds (x, y, z): l_box[i], h_box[i], then "xlo xhi"/"ylo yhi"/"zlo zhi"
    l_box = [0.0]*3
    h_box = [0.0]*3
    for i in range(3):
        l_box[i] = _next(tokens, float)
        h_box[i] = _next(tokens, float)
        _ = _next(tokens)  # xlo/ylo/zlo
        _ = _next(tokens)  # xhi/yhi/zhi

    # Likely "Masses"
    _ = _next(tokens)

    # Mass table from input (indexed by input atom type - 1)
    mass0 = [0.0]*max(n_type0, 1)
    for i in range(n_type0):
        _ = _next(tokens)
        mass0[i] = _next(tokens, float)

    # Likely "Atoms # full"
    _ = _next(tokens)
    _ = _next(tokens)
    _ = _next(tokens)

    # Prepare arrays for Atoms section (only what we need to write)
    id0 = [0]*kkk
    mol_id0 = [0]*kkk
    atom_type0 = [0]*kkk
    q0 = [0.0]*kkk
    x0 = [[0.0, 0.0, 0.0] for _ in range(kkk)]

    # Read kkk atom lines:
    # id mol atom_type q x y z pp pp pp
    # The trailing three floats (image flags etc.) are read and discarded.
    for i in range(kkk):
        id0[i]       = _next(tokens, int)
        mol_id0[i]   = _next(tokens, int)
        atom_type0[i]= _next(tokens, int)
        q0[i]        = _next(tokens, float)
        x0[i][0]     = _next(tokens, float)
        x0[i][1]     = _next(tokens, float)
        x0[i][2]     = _next(tokens, float)
        _ = _next(tokens, float)  # pp
        _ = _next(tokens, float)  # pp
        _ = _next(tokens, float)  # pp

    try:
        with open(OUTPUT_FILE, "w", encoding="utf-8") as out:
            out.write("LAMMPS data file via write_data, version 10 Aug 2015, timestep = 0\n")
            out.write("\n")
            out.write(f"{kkk} atoms\n")
            out.write(f"8 atom types\n")
            out.write("\n")
            # 4 spaces between bounds, then labels
            out.write(f"{l_box[0]}    {h_box[0]} xlo xhi\n")
            out.write(f"{l_box[1]}    {h_box[1]} ylo yhi\n")
            out.write(f"{l_box[2]}    {h_box[2]} zlo zhi\n")
            out.write("\n")
            out.write("Masses\n")
            out.write("\n")
            # Hard-coded masses block
            out.write("1  12.01\n")
            out.write("2   1.008\n")
            out.write("3   14.01\n")
            out.write("4   16\n")
            out.write("5   32.06\n")
            out.write("6   30.97\n")
            out.write("7   35.45\n")
            out.write("8   19\n")
            out.write("\n")
            out.write("Atoms # full\n")
            out.write("\n")

            # Mapping mass -> new atom type
            def map_mass_to_type(m):
                if m == 12.01:
                    return 1
                elif m == 1.008:
                    return 2
                elif m == 14.01:
                    return 3
                elif m == 16:
                    return 4
                elif m == 32.06:
                    return 5
                elif m == 30.97:
                    return 6
                elif m == 35.45:
                    return 7
                elif m == 19:
                    return 8
                else:
                    return None

            for i in range(kkk):
                atype_in = atom_type0[i]
                m = mass0[atype_in - 1] if 1 <= atype_in <= len(mass0) else None
                new_type = map_mass_to_type(m) if m is not None else None

                if new_type is None:
                    print(f"Warning: Unknown mass for atom {id0[i]} - defaulting to type 1", file=sys.stderr)
                    new_type = 1

                out.write(
                    f"{id0[i]}   1   {new_type}   {q0[i]}   {x0[i][0]}   {x0[i][1]}   {x0[i][2]} 0  0  0 \n"
                )
    except Exception as e:
        print(f"system.data_modified could not be opened: {e}", file=sys.stderr)
        sys.exit(1)

    print("File processing completed successfully.")

if __name__ == "__main__":
    main()

