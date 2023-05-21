genome="hg38"

wget -O ${genome}.tar "https://storage.googleapis.com/btc-dshub-pipelines/scRNA/refData/GRCh38.tar"
mkdir ./indexes
tar -xvf ${genome}.tar -C ./indexes
rm -rf ${genome}.tar