#!/bin/sh

# Check what files are actually needed to compile main.c (containing dgesv)

if ! gcc -I include/ main.c blas/*.c lapack/*.c libf2c/*.c -lm; then
  exit 1
fi

for f in blas/*.c lapack/*.c libf2c/*.c; do
  rm "$f"
  if ! gcc -I include/ main.c blas/*.c lapack/*.c libf2c/*.c -lm >& /dev/null; then
    git checkout "$f"
    echo "Need $f"
  else
    echo "Can remove $f"
  fi
done
