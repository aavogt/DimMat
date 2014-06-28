#!bin/bash
set -e
if [[ ! (-d "HList-0.3.4.1") ]]; then
    cabal unpack HList
fi
if [[ ! (-d "numtype-1.1") ]]; then
    cabal unpack numtype
fi
if [[ ! (-d "dimensional-0.13") ]]; then
    cabal unpack dimensional
fi
# https://github.com/feuerbach/standalone-haddock
standalone-haddock --hyperlink-source --package-db=$HOME/.ghc/x86_64-linux-7.8.2/package.conf.d -o doc . HList-* numtype-* dimensional-*
rsync -av --delete-after doc/ ../aavogt.github.io/haddock/DimMat
cd ../aavogt.github.io
git add -A haddock
git commit -m "DimMat $(date)"
git push --mirror
