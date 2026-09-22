#!/bin/bash
# Descarga robusta del snapshot JHTDB con reintentos y resume.
F=/home/delorien/sddf/datos/jhtdb_cache/isotropic1024-coarse-velocity.h5
URL="https://huggingface.co/datasets/ArielLubonja/johns-hopkins-turbulence-database/resolve/main/isotropic1024-coarse-velocity.h5"
mkdir -p "$(dirname "$F")"
for i in 1 2 3 4 5 6 7 8 9 10; do
  echo "[try $i] $(date -Is)"
  curl -C - -L --connect-timeout 30 --speed-time 60 --speed-limit 10000 \
       -o "$F" "$URL" && break
  sleep 20
done
ls -la "$F"