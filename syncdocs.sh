#!bin/bash
set -e
if [[ ! (-d "HList-0.3.0.1") ]]; then
    cabal unpack HList
fi
if [[ ! (-d "numtype-tf-0.1.1") ]]; then
    cabal unpack numtype-tf
fi
if [[ ! (-d "dimensional-tf-0.2.1") ]]; then
    cabal unpack dimensional-tf
fi
# https://github.com/feuerbach/standalone-haddock
standalone-haddock --hyperlink-source --package-db=$HOME/.ghc/x86_64-linux-7.6.2/package.conf.d -o doc . HList-* numtype-tf-* dimensional-tf-*
rsync -av --delete-after doc/ ../aavogt.github.io/haddock/DimMat
cd ../aavogt.github.io
git add -A haddock
git commit -m 'DimMat $(date)'
git push --mirror
