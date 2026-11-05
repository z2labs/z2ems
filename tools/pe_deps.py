#!/usr/bin/env python3
""" Minimal PE import walker (no pefile): lists the DLLs a Windows executable needs, resolved against the given
    folders, recursively.  python pe_deps.py <exe> <dir1> [dir2 ...]   -> prints "<name> <resolved path>" lines,
    and "<name> SYSTEM" for DLLs found nowhere (assumed to be Windows system libraries). """
import os, struct, sys


def imports(path):
    with open(path, 'rb') as fh:
        d = fh.read()
    pe = struct.unpack_from('<I', d, 0x3c)[0]
    nsec = struct.unpack_from('<H', d, pe + 6)[0]
    opt = pe + 24
    magic = struct.unpack_from('<H', d, opt)[0]
    dd = opt + (112 if magic == 0x20b else 96)
    rva, size = struct.unpack_from('<II', d, dd + 8)   # import table = data directory 1
    secs = []
    soff = opt + struct.unpack_from('<H', d, pe + 20)[0]
    for i in range(nsec):
        s = soff + 40 * i
        va, vsz, raw, rsz = struct.unpack_from('<IIII', d, s + 12)[0], struct.unpack_from('<I', d, s + 8)[0], \
            struct.unpack_from('<I', d, s + 20)[0], struct.unpack_from('<I', d, s + 16)[0]
        secs.append((va, max(vsz, rsz), raw))

    def r2o(r):
        for va, sz, raw in secs:
            if va <= r < va + sz:
                return r - va + raw
        return None

    out = []
    if not rva:
        return out
    o = r2o(rva)
    while o is not None:
        name_rva = struct.unpack_from('<I', d, o + 12)[0]
        if name_rva == 0:
            break
        no = r2o(name_rva)
        name = d[no:d.index(b'\0', no)].decode('ascii', 'replace')
        out.append(name)
        o += 20
    return out


def main():
    exe, dirs = sys.argv[1], sys.argv[2:]
    seen, todo, res = {}, [exe], {}
    while todo:
        p = todo.pop()
        for n in imports(p):
            k = n.lower()
            if k in seen:
                continue
            found = None
            for dd in dirs:
                c = os.path.join(dd, n)
                if os.path.exists(c):
                    found = c
                    break
            seen[k] = found or 'SYSTEM'
            if found:
                todo.append(found)
    for k in sorted(seen):
        print(k, seen[k])


if __name__ == '__main__':
    main()
