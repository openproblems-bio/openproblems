export R_LIBS_USER="$HOME/R/Library"
echo ".libPaths(c('$R_LIBS_USER', .libPaths()))" >> $HOME/.Rprofile
pip install --no-cache-dir --user -q -U $1
