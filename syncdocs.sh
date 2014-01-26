#!bin/bash
# https://github.com/feuerbach/standalone-haddock
standalone-haddock --package-db=$HOME/.ghc/x86_64-linux-7.6.2/package.conf.d -o doc .
rsync -av --delete-after doc/ ../aavogt.github.io/haddock/DimMat
cd ../aavogt.github.io
git add -A haddock
git commit -m 'DimMat $(date)'
git push --mirror
