echo "============================================"
echo "Downloading and Installing Miniconda for"
echo "============================================"

export MINICONDA_DOWNLOAD=$MINICONDA_PATH/download
echo "MINICONDA_DOWNLOAD = $MINICONDA_DOWNLOAD"

mkdir -p $MINICONDA_DOWNLOAD;
echo "downloading miniconda.sh for osx ==========="
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O $MINICONDA_DOWNLOAD/miniconda.sh;

echo "installing miniconda ======================="
bash $MINICONDA_DOWNLOAD/miniconda.sh -b -u -p $MINICONDA_PATH;

echo "============================================"
echo "Finished Installing Miniconda"
echo "============================================"
