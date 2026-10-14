"""Compare digit strings: compare.py a.json b.json|bfile key ..."""
import json, sys
def digs(src, key):
    if src.endswith(".txt"):
        return "".join(l.split()[1] for l in open(src) if l.strip() and not l.startswith("#"))
    s = json.load(open(src))[key]
    return s.replace(".", "").lstrip("-")
def agree(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y: break
        n += 1
    return n
if __name__ == "__main__":
    a, b = sys.argv[1], sys.argv[2]
    for key in sys.argv[3:]:
        A, B = digs(a, key), digs(b, key)
        print(f"{key}: agree {agree(A, B)} digits (lens {len(A)}, {len(B)})")
