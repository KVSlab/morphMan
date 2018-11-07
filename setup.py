import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="morphman",
    version="0.1",
    license="GPL",
    author="Aslak W. Bergersen, Henrik A. Kjeldsberg",
    author_email="henrik.kjeldsberg@live.no",
    description="morphman - morphlogical manipulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KVSlab/morphman.git",
    packages=setuptools.find_packages(exclude=["test"]),
    project_urls={
        "Documentation": "https://morphman.readthedocs.io/",
        "Source Code": "https://github.com/KVSlab/morphman",
    },
    install_requires=[
        'scipy', 'numpy', 'vtk'
    ],
    dependency_links=[
        'https://github.com/vmtk/vmtk'
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        "Programming Language :: Python :: 3",
    ],
)
