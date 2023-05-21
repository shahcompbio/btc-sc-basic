# docker build --platform linux/amd64 -t scaligners:2.0 -f container/scAligners.Dockerfile .
docker build --platform linux/amd64 -t scpackages:2.0 -f container/scPackages.Dockerfile .

docker tag scpackages:2.0 zatzmanm/scpackages:2.0
docker push zatzmanm/scpackages:2.0