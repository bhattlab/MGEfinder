cp ../env/conda_linux64.yaml environment.yml
docker build -t mgefinder .
rm environment.yml
