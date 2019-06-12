rm -r mgefinder.egg-info/ dist/ build/
python setup.py sdist
python setup.py bdist_wheel --universal
twine upload dist/*
