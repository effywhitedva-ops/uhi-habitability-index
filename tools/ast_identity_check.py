#!/usr/bin/env python3
"""Check that commenting changed no code: compare ast.dump of each module with docstrings stripped.

Usage: python tools/ast_identity_check.py <original_root> <commented_root> file1.py file2.py ...
Exit code 0 only if every listed module is AST-identical once docstrings are removed.
"""
import ast, sys, os

def strip_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            b = node.body
            if b and isinstance(b[0], ast.Expr) and isinstance(b[0].value, ast.Constant) and isinstance(b[0].value.value, str):
                node.body = b[1:] or [ast.Pass()]
    return tree

def dump(path):
    return ast.dump(strip_docstrings(ast.parse(open(path, encoding="utf-8").read())), include_attributes=False)

orig, new, files = sys.argv[1], sys.argv[2], sys.argv[3:]
bad = 0
for f in files:
    same = dump(os.path.join(orig, f)) == dump(os.path.join(new, f))
    print(("IDENTICAL " if same else "DIFFERENT ") + f); bad += not same
print("ALL_IDENTICAL" if bad == 0 else f"{bad} MODULE(S) DIFFER"); sys.exit(1 if bad else 0)
