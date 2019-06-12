from setuptools import setup, find_packages

setup(
    name="mgefinder",
    version='0.0.1.dev7',
    description='A toolbox for identifying mobile genetic element (MGE) insertions from short-read sequencing data of bacterial isolates.',
    url='https://github.com/bhattlab/MGEfinder',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click',
        'pandas',
        'biopython',
        'pysam',
        'editdistance',
        'scipy',
        'networkx',
        'tqdm'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'mgefinder = mgefinder.main:cli'
        ],
}
)
