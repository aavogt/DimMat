#!bin/bash
set -e
if [[ ! (-d "HList-0.3.0.1") ]]; then
    cabal unpack HList
fi
# https://github.com/feuerbach/standalone-haddock
standalone-haddock --package-db=$HOME/.ghc/x86_64-linux-7.6.2/package.conf.d -o doc . HList-*
rsync -av --delete-after doc/ ../aavogt.github.io/haddock/DimMat
cd ../aavogt.github.io
git add -A haddock
git commit -m 'DimMat $(date)'
git push --mirror
