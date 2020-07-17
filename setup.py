import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

DEPENDENCIES = []  # 'scipy', 'numpy', 'vtk']
TEST_DEPENDENCIES = []  # 'pytest']

VERSION = "1.0"
URL = "https://github.com/KVSlab/morphMan.git"

setuptools.setup(
    name="morphman",
    version=VERSION,
    license="GPL",
    author="Aslak W. Bergersen, Henrik A. Kjeldsberg",
    author_email="henrik.kjeldsberg@live.no",
    url=URL,
    project_urls={
        "Documentation": "https://morphman.readthedocs.io/",
        "Source Code": URL,
    },
    description="morphman - morphlogical manipulation",
    long_description=long_description,
    long_description_content_type="text/markdown",

    # Dependencies
    install_requires=DEPENDENCIES,
    tests_require=TEST_DEPENDENCIES,
    dependency_links=['https://github.com/vmtk/vmtk'],

    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        "Programming Language :: Python :: 3",
    ],
    packages=["morphman",
              "morphman.common",
              "morphman.misc",
              "morphman.automated_landmarking"],
    package_dir={"morphman": "morphman"},
    entry_points={'console_scripts':
                      ['morphman-area=morphman.manipulate_area:main_area',
                       'morphman-bend=morphman.manipulate_bend:main_bend',
                       'morphman-bifurcation=morphman.manipulate_bifurcation:main_bifurcation',
                       'morphman-curvature=morphman.manipulate_curvature:main_curvature',
                       'morphman-branch=morphman.manipulate_branch:main_branch',
                       'morphman-surface=morphman.manipulate_surface:main_surface']
                  }
)
