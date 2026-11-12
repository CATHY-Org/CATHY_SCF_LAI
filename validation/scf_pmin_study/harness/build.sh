#!/usr/bin/env bash
# Build the two CATHY executables used in the SCF/PMIN smoke tests.
#   cathy_classic : CATHY-Org/CATHY_src_org @ pinned commit (SCF = one global scalar)
#   cathy_scf     : CATHY-Org/CATHY_SCF_LAI main        (SCF per vegetation type + joint cap)
# Both are built from *scratch copies*; upstream sources are never modified except for:
#   (1) a symlink conf_fac.f -> Conf_Fac.f (Makefile/filename case mismatch on Linux)
#   (2) silencing per-node debug WRITE(6,*) in the SCF build's etran.f (44 MB stdout/hour otherwise)
set -euo pipefail
ROOT="${1:-$(pwd)}"; BIN="$ROOT/bin"; mkdir -p "$BIN" "$ROOT/_build"
CLASSIC_URL=https://github.com/CATHY-Org/CATHY_src_org.git
CLASSIC_COMMIT=ef54607df776dbb6f7d974329e5f72a276403d5c
SCF_URL=https://github.com/CATHY-Org/CATHY_SCF_LAI.git

build () {  # $1=name  $2=src dir
  local d="$ROOT/_build/$1"; rm -rf "$d"; cp -r "$2" "$d"; cd "$d"
  rm -f *.o cathy_ft
  [ -e Conf_Fac.f ] && ln -sf Conf_Fac.f conf_fac.f
  if [ "$1" = scf ]; then python3 - <<'PY'
import re
s=open('etran.f').read().split('\n'); out=[]; skip=0
for i,l in enumerate(s):
    if skip: skip-=1; continue
    if re.match(r"\s*write\(6,\*\)\s*'(I=|ZROOT=)",l, re.I):
        # 'I=' BTRAN debug write has one continuation line
        if "BTRAN=" in l: skip=1
        continue
    out.append(l)
open('etran.f','w').write('\n'.join(out))
PY
  fi
  make -j8 >/dev/null 2>"$ROOT/_build/$1.log" || { tail -20 "$ROOT/_build/$1.log"; exit 1; }
  cp cathy_ft "$BIN/cathy_$1"
}

if [ ! -d "$ROOT/_src/classic" ]; then
  mkdir -p "$ROOT/_src/classic"; cd "$ROOT/_src/classic"; git init -q
  git remote add origin $CLASSIC_URL; git fetch -q --depth 1 origin $CLASSIC_COMMIT; git checkout -q FETCH_HEAD
fi
[ -d "$ROOT/_src/scf" ] || git clone -q --depth 1 $SCF_URL "$ROOT/_src/scf"
build classic "$ROOT/_src/classic/src"
build scf     "$ROOT/_src/scf/src"
( cd "$ROOT/_src/scf" && echo "scf_commit=$(git rev-parse --short HEAD)" ) > "$BIN/VERSIONS.txt"
echo "classic_commit=$CLASSIC_COMMIT" >> "$BIN/VERSIONS.txt"
ls -la "$BIN"
