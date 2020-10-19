from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()

setup(
    name="mgefinder",
    version='1.0.6',
    description='A toolbox for identifying mobile genetic element (MGE) insertions from short-read sequencing data of bacterial isolates.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/bhattlab/MGEfinder',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click==7.0',
        'pandas==0.25.3',
        'biopython==1.75',
        'pysam==0.15.3',
        'editdistance',
        'scipy==1.4.0',
        'networkx==2.4',
        'tqdm==4.40.2'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'mgefinder = mgefinder.main:cli'
        ]
}
)
