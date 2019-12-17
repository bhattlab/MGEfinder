cp ../env/conda_linux64.yaml .
docker build -t mgefinder .
rm conda_linux64.yaml
